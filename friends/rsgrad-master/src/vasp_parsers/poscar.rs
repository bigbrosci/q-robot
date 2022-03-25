use std::{
    cmp::Ordering,
    path::Path,
    fs,
    io::{
        BufReader, 
        BufRead,
    },
    fmt,
};
use anyhow::{anyhow, Context};
use crate::{
    Result,
    Structure,
    Mat33,
    MatX3,
};

#[derive(Clone, Debug)]
pub struct Poscar {  // I have no plan to support vasp4 format
    pub comment: String,
    pub scale: f64,
    pub cell: Mat33<f64>,
    pub ion_types: Vec<String>,
    pub ions_per_type: Vec<i32>,
    pub pos_cart: MatX3<f64>,
    pub pos_frac: MatX3<f64>,
    pub constraints: Option<MatX3<bool>>,
}


impl Poscar {
    pub fn from_file(path: &(impl AsRef<Path> + ?Sized)) -> Result<Self> {
        //  Read to the first emtpy line then parse it.
        let f = fs::File::open(path)?;
        let txt = BufReader::new(f).lines()
            .take_while(|x| x.is_ok())
            .map(|x| x.unwrap())
            .take_while(|x| !x.trim().is_empty())
            .collect::<Vec<_>>()
            .join("\n");
        Self::from_str(&txt)
    }

    pub fn from_str(txt: &str) -> Result<Self> {
        let mut lines = txt.lines();
        let comment: String = lines.next().context("[POSCAR]: File may be blank.")?.trim().to_string();
        let scale: f64 = lines.next().context("[POSCAR]: Cannot parse scale constant.")?
            .split_whitespace()
            .next().context("[POSCAR]: Scale line may be empty.")?
            .parse::<f64>()
            .context("[POSCAR]: Scale constant cannot be converted to float number.")?;
        
        let scale = match scale.partial_cmp(&0.0) {
            Some(Ordering::Greater) => scale,
            Some(Ordering::Less) | Some(Ordering::Equal) => 
                return Err(anyhow!("[POSCAR]: Scale constant should be greater than 0.0.")),
            None => return Err(anyhow!("[POSCAR]: Scale constant cannot be NaN.")),
        };
        
        let cell: Mat33<f64> = {
            let mut v = [[0.0f64; 3]; 3];
            for i in 0..3 {
                let line = lines.next().context("[POSCAR]: Incomplete lines for cell info.")?;
                let row = line.split_whitespace().take(3).collect::<Vec<_>>();
                if row.len() < 3 {
                    return Err(anyhow!("[POSCAR]: Cell lines incomplete."));
                }
                for (j, x) in row.into_iter().enumerate() {
                    let val = x.parse::<f64>().context(format!("[POSCAR]: Cell lines contain invalid value: `{}` .", x))?;
                    if val.is_nan() {
                        return Err(anyhow!("[POSCAR]: Cell lines contain NaN value."));
                    }
                    v[i][j] = val;
                }
            }
            v
        };
        
        let ion_types = {
            let words = lines.next()
                .context("[POSCAR]: Element tags line not found, rsgrad has no plan to support vasp4 format.")?
                .split_whitespace()
                .take_while(|x| !x.contains("!"))
                .map(|x| x.to_string())
                .collect::<Vec<_>>();
            if words.len() == 0 {
                return Err(anyhow!("[POSCAR]: At lease one element is needed."));
            }
            words
        };

        let ions_per_type = {
            let numbers = lines.next()
                .context("[POSCAR]: Count of each element not found.")?
                .split_whitespace()
                .take_while(|x| !x.contains("!"))
                .map(|x| x.parse::<i32>().context("[POSCAR]: Invalid atom count of element."))
                .collect::<Result<Vec<_>>>()?;
            if numbers.len() != ion_types.len() {
                return Err(anyhow!("[POSCAR]: Inconsistent element types and atom counts."));
            }
            if numbers.iter().any(|x| x <= &0) {
                return Err(anyhow!("[POSCAR]: Atom counts cannot be zero or negative."));
            }
            numbers
        };

        let line = lines.next().context("[POSCAR]: Constraints or Coordination type not found.")?;
        let has_constraints = {
            match line.trim_start().chars().next() {
                Some('s') | Some('S') => true,
                Some('d') | Some('D') | Some('c') | Some('C') | Some('k') | Some('K') => false,
                _ => return Err(anyhow!("[POSCAR]: Constraints line or Coordination type line missing."))
            }
        };
        
        let is_direct = {
            let line = if has_constraints {
                lines.next().context("[POSCAR]: Coordination type not found.")?
            } else {
                line
            };
            match line.trim_start().chars().next() {
                Some('d') | Some('D') => true,
                Some('c') | Some('C') | Some('k') | Some('K') => false,
                _ => return Err(anyhow!("[POSCAR]: Coordination type line missing."))
            }
        };

        let mut coords: MatX3<f64> = vec![];
        let mut constraints: Option<MatX3<bool>> = if has_constraints { Some(vec![]) } else { None };
        while let Some(line) = lines.next() {
            if line.trim().is_empty() {
                break;
            }
            let v = line.split_whitespace().collect::<Vec<_>>();
            coords.push( [ v[0].parse::<f64>().context(format!("[POSCAR]: Coordinates value invalid: `{}` .", v[0]))?,
                           v[1].parse::<f64>().context(format!("[POSCAR]: Coordinates value invalid: `{}` .", v[1]))?,
                           v[2].parse::<f64>().context(format!("[POSCAR]: Coordinates value invalid: `{}` .", v[2]))?, ]);

            if let Some(c) = &mut constraints {
                let mut _c = [true; 3];
                for i in 0..3 {
                    _c[i] = match v[i+3].chars().next().unwrap() {
                        't' | 'T' => true,
                        'f' | 'F' => false,
                        _ => return Err(anyhow!("[POSCAR]: Constraints should be either 'T' or 'F'")),
                    };
                }
                c.push(_c);
            }
        }

        if coords.len() as i32 != ions_per_type.iter().sum::<i32>() {
            return Err(anyhow!("[POSCAR]: Count of coordinates inconsistent with sum of atom counts."));
        }

        let (pos_cart, pos_frac): (MatX3<f64>, MatX3<f64>) = {
            if is_direct {
                let cart = Self::convert_frac_to_cart(&coords, &cell);
                (cart, coords)
            } else {
                let frac = Self::convert_cart_to_frac(&coords, &cell)
                    .context("[POSCAR]: Cell matrix is singular.")?;
                (coords, frac)
            }
        };

        // TODO parse velocity, may be implemented later, if needed.

        Ok(Poscar{
            comment,
            scale,
            cell,
            ion_types,
            ions_per_type,
            pos_cart,
            pos_frac,
            constraints
        })
    }


    pub fn from_structure(s: Structure) -> Self {
        Self {
            comment: "Generated by rsgrad".to_string(),
            scale: 1.0,
            cell: s.cell,
            ion_types: s.ion_types,
            ions_per_type: s.ions_per_type,
            pos_cart: s.car_pos,
            pos_frac: s.frac_pos,
            constraints: s.constr,
        }
    }


    pub fn into_structure(self) -> Structure {
        let self2 = self.normalize();

        Structure {
            cell: self2.cell,
            ion_types: self2.ion_types,
            ions_per_type: self2.ions_per_type,
            car_pos: self2.pos_cart,
            frac_pos: self2.pos_frac,
            constr: self2.constraints,
        }
    }


    pub fn to_formatter(&self) -> PoscarFormatter<'_> {
        PoscarFormatter::new(&self)
    }


    pub fn get_cell_params(&self) -> ([f64; 3], [f64; 3]) {
        let lengths = {
            let a = &self.cell;
            let mut ret = [0.0; 3];
            ret[0] = f64::sqrt(a[0][0] * a[0][0] + a[0][1] * a[0][1] + a[0][2] * a[0][2]) * self.scale;
            ret[1] = f64::sqrt(a[1][0] * a[1][0] + a[1][1] * a[1][1] + a[1][2] * a[1][2]) * self.scale;
            ret[2] = f64::sqrt(a[2][0] * a[2][0] + a[2][1] * a[2][1] + a[2][2] * a[2][2]) * self.scale;
            ret
        };

        let product = |a: &[f64; 3], b: &[f64; 3]| {
            a[0] * b[0]
          + a[1] * b[1]
          + a[2] * b[2]
        };

        let angles = {
            let a = &self.cell;
            let alpha = f64::acos(product(&a[1], &a[2]) * self.scale.powi(2) / (lengths[1] * lengths[2])) * 
                        180.0 / std::f64::consts::PI;
            let beta  = f64::acos(product(&a[0], &a[2]) * self.scale.powi(2) / (lengths[0] * lengths[2])) * 
                        180.0 / std::f64::consts::PI;
            let gamma = f64::acos(product(&a[0], &a[1]) * self.scale.powi(2) / (lengths[0] * lengths[1])) * 
                        180.0 / std::f64::consts::PI;
            [alpha, beta, gamma]
        };

        (lengths, angles)
    }


    pub fn normalize(mut self) -> Self {
        for i in 0..3 {
            for j in 0..3 {
                self.cell[i][j] *= self.scale;
            }
        }

        for i in 0..self.pos_cart.len() {
            for j in 0..3 {
                self.pos_cart[i][j] *= self.scale;
            }
        }
        self.scale = 1.0;

        self
    }


    pub fn get_natoms(&self) -> i32 {
        self.ions_per_type.iter().sum()
    }


    pub fn get_ntypes(&self) -> i32 {
        self.ion_types.len() as i32
    }


    pub fn get_volume(&self) -> f64 {
        Self::mat33_det(&self.cell) * self.scale.powi(3)
    }


    #[inline]
    fn mat33_det(mat: &Mat33<f64>) -> f64 {
        return mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2])
             - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])
             + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
    }


    fn mat33_inv(mat: &Mat33<f64>) -> Option<Mat33<f64>> {
        let det = Self::mat33_det(mat);
        if det.abs() < 1E-5 { return None; }
        let invdet = 1.0 / det;
        let mut bmat: Mat33<f64> = [[0.0; 3]; 3];
        bmat[0][0] = (mat[1][1] * mat[2][2]  -  mat[2][1] * mat[1][2]) * invdet;
        bmat[0][1] = (mat[0][2] * mat[2][1]  -  mat[0][1] * mat[2][2]) * invdet;
        bmat[0][2] = (mat[0][1] * mat[1][2]  -  mat[0][2] * mat[1][1]) * invdet;
        bmat[1][0] = (mat[1][2] * mat[2][0]  -  mat[1][0] * mat[2][2]) * invdet;
        bmat[1][1] = (mat[0][0] * mat[2][2]  -  mat[0][2] * mat[2][0]) * invdet;
        bmat[1][2] = (mat[1][0] * mat[0][2]  -  mat[0][0] * mat[1][2]) * invdet;
        bmat[2][0] = (mat[1][0] * mat[2][1]  -  mat[2][0] * mat[1][1]) * invdet;
        bmat[2][1] = (mat[2][0] * mat[0][1]  -  mat[0][0] * mat[2][1]) * invdet;
        bmat[2][2] = (mat[0][0] * mat[1][1]  -  mat[1][0] * mat[0][1]) * invdet;
        Some(bmat)
    }


    fn matx3_mul_mat33(matx3: &MatX3<f64>, mat33: &Mat33<f64>) -> MatX3<f64> {
        let len = matx3.len();
        let mut ret = vec![[0.0; 3]; len];
        for i in 0..len {
            // manual loop unroll
            ret[i][0] += matx3[i][0] * mat33[0][0];
            ret[i][0] += matx3[i][1] * mat33[1][0];
            ret[i][0] += matx3[i][2] * mat33[2][0];

            ret[i][1] += matx3[i][0] * mat33[0][1];
            ret[i][1] += matx3[i][1] * mat33[1][1];
            ret[i][1] += matx3[i][2] * mat33[2][1];

            ret[i][2] += matx3[i][0] * mat33[0][2];
            ret[i][2] += matx3[i][1] * mat33[1][2];
            ret[i][2] += matx3[i][2] * mat33[2][2];
        }
        ret
    }


    pub fn convert_cart_to_frac(cart: &MatX3<f64>, cell: &Mat33<f64>) -> Option<MatX3<f64>> {
        let inv = Self::mat33_inv(cell)?;
        Some(Self::matx3_mul_mat33(cart, &inv))
    }


    pub fn convert_frac_to_cart(frac: &MatX3<f64>, cell: &Mat33<f64>) -> MatX3<f64> {
        Self::matx3_mul_mat33(frac, cell)
    }
}


impl From<Structure> for Poscar {
    fn from(s: Structure) -> Self {
        Self::from_structure(s)
    }
}


pub struct PoscarFormatter<'a> {
    pub poscar: &'a Poscar,
    pub preserve_constraints: bool,
    pub fraction_coordinates: bool,
    pub add_symbol_tags: bool,
}


impl<'a> PoscarFormatter<'a> {
    pub fn new(poscar: &'a Poscar) -> Self {
        Self {
            poscar: &poscar,
            preserve_constraints: true,
            fraction_coordinates: true,
            add_symbol_tags: true,
        }
    }

    pub fn preserve_constraints(mut self, flag: bool) -> Self {
        self.preserve_constraints = flag;
        self
    }

    pub fn fraction_coordinates(mut self, flag: bool) -> Self {
        self.fraction_coordinates = flag;
        self
    }

    pub fn add_symbol_tags(mut self, flag: bool) -> Self {
        self.add_symbol_tags = flag;
        self
    }

    pub fn to_file(&self, path: &(impl AsRef<Path> + ?Sized)) -> Result<()> {
        fs::write(path, self.to_string())?;
        Ok(())
    }
}


impl<'a> fmt::Display for PoscarFormatter<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let poscar = &self.poscar;

        writeln!(f, "{}", &poscar.comment)?;
        writeln!(f, "{:10.7}", poscar.scale)?;

        for i in 0..3 {
            writeln!(f, "   {:15.9}   {:15.9}   {:15.9}", poscar.cell[i][0], poscar.cell[i][1], poscar.cell[i][2])?;
        }

        {
            let mut symbol_line = String::with_capacity(8);
            let mut count_line = String::with_capacity(8);

            for (t, c) in poscar.ion_types.iter().zip(poscar.ions_per_type.iter()) {
                symbol_line += &format!(" {:>6}", t);
                count_line += &format!(" {:>6}", c);
            }

            write!(f, "{}\n{}\n", symbol_line, count_line)?;
        }

        let atom_symbol_index = {
            let mut ret = vec![String::new(); 0];
            let mut ind = 0;
            
            for (symbol, count) in poscar.ion_types.iter().zip(poscar.ions_per_type.iter()) {
                for i in 1..=*count {
                    ind += 1;
                    ret.push(format!("{:>6}-{:03}  {:3}", symbol, i, ind));
                }
            }
            ret
        };

        let (write_constraints, constr) = if poscar.constraints.is_some() && self.preserve_constraints {
            f.write_str("Selective Dynamics\n")?;
            (true, poscar.constraints.as_ref().unwrap().clone())
        } else {
            (false, vec![[true, true, true]; 0])
        };

        let coords = if self.fraction_coordinates {
            f.write_str("Direct\n")?;
            &poscar.pos_frac
        } else {
            f.write_str("Cartesian\n")?;
            &poscar.pos_cart
        };

        for i in 0..coords.len() {
            write!(f, "  {:16.10}  {:16.10}  {:16.10} ", coords[i][0], coords[i][1], coords[i][2])?;

            if write_constraints {
                for c in constr[i] {
                    f.write_str(if c { "  T " } else { "  F " })?;
                }
            }

            if self.add_symbol_tags {
                write!(f, "! {}", &atom_symbol_index[i])?;
            }
            writeln!(f)?;
        }

        Ok(())
    }
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mat33_det() {
        let mat = [[ 1.0, 2.0, 3.0], 
                   [ 4.0, 5.0, 6.0],
                   [ 7.0, 8.0, 9.0]];
        assert_eq!(Poscar::mat33_det(&mat), 0.0);

        let mat = [[ 1.0, 3.0, 4.0], 
                   [ 2.0, 1.0, 5.0],
                   [ 5.0, 3.0, 0.0]];
        assert_eq!(Poscar::mat33_det(&mat), 64.0);
    }

    #[test]
    fn test_mat33_inv() {
        let mat = [[ 1.0, 2.0, 3.0], 
                   [ 4.0, 5.0, 6.0],
                   [ 7.0, 8.0, 9.0]];
        assert_eq!(Poscar::mat33_inv(&mat), None);

        let mat = [[ 1.0, 3.0, 4.0], 
                   [ 2.0, 1.0, 5.0],
                   [ 5.0, 3.0, 0.0]];
        assert_eq!(Poscar::mat33_inv(&mat), 
                   Some([[-0.234375,  0.1875,  0.171875],
                         [ 0.390625, -0.3125,  0.046875],
                         [ 0.015625,  0.1875, -0.078125]]));
    }

    #[test]
    fn test_matx3_mul_mat33() {
        let matx3 = vec![[ 1.0, 2.0, 3.0], 
                         [ 4.0, 5.0, 6.0],
                         [ 7.0, 8.0, 9.0]];
        let mat33 = [[ 1.0, 3.0, 4.0], 
                     [ 2.0, 1.0, 5.0],
                     [ 5.0, 3.0, 0.0]];
        assert_eq!(Poscar::matx3_mul_mat33(&matx3, &mat33), 
                   vec![[20., 14., 14.],
                        [44., 35., 41.],
                        [68., 56., 68.]]);
    }


    #[test]
    fn test_convert_cart_to_frac() {
        let cell: Mat33<f64> = [[7.50591018156692, -4.33353926384082, -1.4e-15], 
                    [0.0, 8.66707852768163, 0.0], 
                    [0.0, 0.0, 66.9999794426203]];
        let frac = vec![[0.9992791851439673, 0.9999905627575514, 0.1300144910859293]];
        let cart = vec![[7.50049981, 4.33658115, 8.71096823]];
        assert_eq!(Poscar::convert_cart_to_frac(&cart, &cell), Some(frac));
    }

    #[test]
    fn test_convert_frac_to_cart() {
        let cell: Mat33<f64> = [[7.50591018156692, -4.33353926384082, -1.4e-15], 
                    [0.0, 8.66707852768163, 0.0], 
                    [0.0, 0.0, 66.9999794426203]];
        let frac = vec![[0.99927918, 0.99999056, 0.13001449]];
        let cart = vec![[7.500499771389843, 4.336581148391669, 8.710968157242762]];
        assert_eq!(Poscar::convert_frac_to_cart(&frac, &cell), cart)
    }
}
