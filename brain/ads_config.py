

def find_rhombus(atoms, top_list, B, D, bond_threshold=2.6):
    """
    Given a list of four atom indices (top_list) that form a quadrilateral in clockwise order,
    this function identifies a valid rhombus if:
      1) The pair of atoms with the longest distance is assigned as A and C (with the smaller index as A).
      2) The remaining two atoms are B and D (with the smaller index as B).
      3) The bonds A–B, B–C, C–D, and D–A all have distances less than or equal to bond_threshold.
         (In this design, B and D are assumed to be bonded, representing the shorter bridge.)
    
    Parameters:
      atoms: An object that supports atoms.get_distance(i, j) for distance between atoms.
      top_list: List of four atom indices.
      B, D: The two atoms chosen (from top_list) that form the short bridge (B–D).
      bond_threshold: Maximum allowed distance for a bond.
    
    Returns:
      A list [A, B, C, D] representing the valid rhombus in clockwise order,
      or an empty list if no valid rhombus is found.
    """
    # Compute the distance between B and D (the short bridge)
    BD_distance = atoms.get_distance(B, D)
    top_list_sorted = sorted(top_list)

    for A in top_list_sorted:
        if A == B or A == D:
            continue

        # Check whether A forms bonds with both B and D
        if atoms.get_distance(A, B) <= bond_threshold and atoms.get_distance(A, D) <= bond_threshold:
            for C in top_list_sorted:
                if C in (A, B, D):
                    continue

                # Check whether C forms bonds with both B and D
                if atoms.get_distance(C, B) <= bond_threshold and atoms.get_distance(C, D) <= bond_threshold:
                    AC_distance = atoms.get_distance(A, C)
                    # The longest distance should be between A and C.
                    if AC_distance > BD_distance:
                        return [A, B, C, D]
    return []


def dihedral_angle(points):
    """
    Calculate the dihedral angle (in degrees) between the plane defined by A, B, D
    and the plane defined by B, C, D.
    
    Parameters:
        points (array-like): A 4x3 array containing the positions of points A, B, C, and D.
    
    Returns:
        float: The dihedral angle in degrees.
    """
    A, B, C, D = points
    # Compute normals for the two planes
    n1 = np.cross(B - A, D - A)
    n1 /= np.linalg.norm(n1)
    n2 = np.cross(C - B, D - B)
    n2 /= np.linalg.norm(n2)
    # Compute the angle between the normals
    dot_val = np.clip(np.dot(n1, n2), -1.0, 1.0)
    angle = np.degrees(np.arccos(dot_val))
    return angle

def filter_rhombus_by_dihedral(rhombuses, atoms, dihedral_threshold=120):
    """
    Separate candidate rhombus configurations into planar and non-planar groups
    based on the dihedral angle between planes ABD and BCD.
    
    Parameters:
        rhombuses (list): List of candidate rhombus configurations (each is a list [A, B, C, D]).
        atoms: An ASE Atoms object (with .position and get_distance() method).
        dihedral_threshold (float): Angle (in degrees) above which a configuration is considered planar.
    
    Returns:
        tuple: (planar, non_planar)
            - planar: list of candidates with dihedral angle >= dihedral_threshold.
            - non_planar: list of candidates with dihedral angle < dihedral_threshold.
    """
    planar = []
    non_planar = []
    for indices in rhombuses:
        points = np.array([atoms[idx].position for idx in indices])
        angle = dihedral_angle(points)
        if angle >= dihedral_threshold:
            planar.append(indices)
        else:
            non_planar.append(indices)
    return planar, non_planar

def filter_rhombus_by_AC_BD(planar, non_planar, atoms, AC_cutoff=3.0):
    """
    From the candidate rhombus configurations (planar and non-planar), check if any candidate in the non-planar group
    has a very short A–C distance (less than AC_cutoff). If so, use the A–C pair (sorted) to find candidates
    that have this pair as their B–D bond. Those candidates are then removed from both the planar and non-planar lists
    and collected in an invalid_rhombus list.
    
    Parameters:
        planar (list): List of candidate rhombus configurations classified as planar.
        non_planar (list): List of candidate configurations classified as non-planar.
        atoms: An ASE Atoms object.
        AC_cutoff (float): Distance cutoff (in Å) for the A–C bond.
    
    Returns:
        tuple: (new_planar, new_non_planar, invalid_rhombus)
            - new_planar: Planar candidates after removal.
            - new_non_planar: Non-planar candidates after removal.
            - invalid_rhombus: Candidates removed because their A–C pair (from non-planar with short A–C)
                                is used as a B–D pair.
    """
    invalid_rhombus = []
    # Combine candidates from both groups
    all_candidates = planar + non_planar

    # Build a dictionary mapping each candidate's B–D pair (sorted) to candidate(s)
    bd_dict = {}
    for candidate in all_candidates:
        BD = tuple(sorted((candidate[1], candidate[3])))
        bd_dict.setdefault(BD, []).append(candidate)

    # For each candidate in the non-planar list, check its A–C distance.
    for candidate in non_planar:
        AC = tuple(sorted((candidate[0], candidate[2])))
        AC_distance = atoms.get_distance(candidate[0], candidate[2])
        if AC_distance < AC_cutoff:
            if AC in bd_dict:
                for cand in bd_dict[AC]:
                    if cand not in invalid_rhombus:
                        invalid_rhombus.append(cand)

    new_planar = [cand for cand in planar if cand not in invalid_rhombus]
    new_non_planar = [cand for cand in non_planar if cand not in invalid_rhombus]

    return new_planar, new_non_planar, invalid_rhombus


def point_in_prism(point, base_vertices, normal, height):
    """
    Check if a point is inside the prism defined by the base triangle and its top translated version.

    Parameters:
    - point (ndarray): The position of the atom being checked.
    - base_vertices (list): The three vertices of the base triangle.
    - normal (ndarray): The normal vector of the triangle.
    - height (float): The prism height limit.

    Returns:
    - bool: True if the point is inside the prism, False otherwise.
    """
    def is_point_in_triangle(p, v1, v2, v3):
        """Check if a point p lies inside a 2D projection of a triangle"""
        u = v2 - v1
        v = v3 - v1
        w = p - v1

        uu = np.dot(u, u)
        uv = np.dot(u, v)
        vv = np.dot(v, v)
        wu = np.dot(w, u)
        wv = np.dot(w, v)

        denom = uv * uv - uu * vv
        if abs(denom) < 1e-6:
            return False  # Degenerate triangle

        s = (uv * wv - vv * wu) / denom
        t = (uv * wu - uu * wv) / denom

        return (s >= 0) and (t >= 0) and (s + t <= 1)

    # Project point onto the base plane
    v0 = base_vertices[0]
    distance_to_base = np.dot(point - v0, normal)

    # **Check if point is between top and bottom planes**
    if not (-height <= distance_to_base <= height):
        return False

    # Project the point onto the base triangle's plane
    projected_point = point - distance_to_base * normal

    # **Check if projected point is inside the base triangle**
    return is_point_in_triangle(projected_point, *base_vertices)

def count_atoms_in_prism(atoms, indices, height):
    """
    Count the number of atoms within the prism formed by translating the base triangle
    along its normal vector by a given height.
    """
    pos = atoms.get_positions()

    # Get positions of the three triangle vertices
    p1, p2, p3 = pos[indices[0]], pos[indices[1]], pos[indices[2]]

    # Compute the normal vector
    normal = calculate_normal(p1, p2, p3)

    # **Translate base vertices to form the top and bottom triangles**
    top_triangle = [p1 + height * normal, p2 + height * normal, p3 + height * normal]
    bottom_triangle = [p1 - height * normal, p2 - height * normal, p3 - height * normal]

    # Count atoms in the top and bottom regions
    atom_count_top = 0
    atom_count_bottom = 0

    for i in range(len(pos)):
        if i in indices:  # Skip the three triangle atoms
            continue

        p = pos[i]

        # **Check if the atom is inside the top prism region**
        if point_in_prism(p, top_triangle, normal, height):
            print(f"Top Prism: {p}")
            atom_count_top += 1

        # **Check if the atom is inside the bottom prism region (now correctly defined)**
        if point_in_prism(p, bottom_triangle, normal, height):  
            atom_count_bottom += 1
            print(f"Bottom Prism: {p}")

    return atom_count_top, atom_count_bottom

def is_exposed_triangle(atoms, triangle_sites, height=2.2):
    """
    Check if a triangle site is exposed based on atom count in the prism region.

    Returns:
        True  -> Triangle is exposed
        False -> Triangle is not exposed (blocked)
    """
    atom_count_top, atom_count_bottom = count_atoms_in_prism(atoms, triangle_sites, height)

    if atom_count_top > 0 and atom_count_bottom > 0:
        # print(atom_count_top, atom_count_bottom)
        return False 
    return True 

def is_exposed_rhombus(atoms, rhombus_indices, height=2.5):
    """
    Check if a rhombus formed by four atom indices is exposed.
    
    The function divides the rhombus (defined by indices [A, B, C, D]) into two triangles:
      - Triangle 1: [A, B, D]
      - Triangle 2: [B, C, D]
    
    For each triangle, it calls count_atoms_in_prism(atoms, triangle, height) to determine the number
    of atoms above (top) and below (bottom) the plane of the triangle. If, for both triangles, at least
    one triangle does NOT have atoms on both sides of its plane, the rhombus is considered exposed.
    
    Parameters:
        atoms: An ASE Atoms object (each atom has a .position attribute).
        rhombus_indices (list): A list of four atom indices [A, B, C, D] (in clockwise order).
        height (float): Height threshold (in Å) for the prism check.
    
    Returns:
        bool: True if the rhombus is exposed, False otherwise.
    """
    triangles = [
        [rhombus_indices[0], rhombus_indices[1], rhombus_indices[3]],  # Triangle ABD
        [rhombus_indices[1], rhombus_indices[2], rhombus_indices[3]]   # Triangle BCD
    ]
    
    results = []
    for triangle in triangles:
        atom_count_top, atom_count_bottom = count_atoms_in_prism(atoms, triangle, height)
        # If both top and bottom have atoms, then the triangle is not exposed.
        if atom_count_top > 0 and atom_count_bottom > 0:
            results.append(False)
        else:
            results.append(True)
    return all(results)


def filter_exposed_rhombus(candidates, atoms, height=2.5):
    """
    Filter candidate rhombus configurations based on their exposure.
    
    For each candidate (a list of four atom indices), this function uses is_exposed_rhombus()
    to determine whether the configuration is exposed.
    
    Parameters:
        candidates (list): List of candidate rhombus configurations (each is a list [A, B, C, D]).
        atoms: An ASE Atoms object.
        height (float): Height threshold (in Å) for exposure.
    
    Returns:
        tuple: (exposed_candidates, not_exposed_candidates)
            - exposed_candidates: List of candidates that are exposed.
            - not_exposed_candidates: List of candidates that are not exposed.
    """
    exposed = []
    not_exposed = []
    for candidate in candidates:
        if is_exposed_rhombus(atoms, candidate, height):
            exposed.append(candidate)
        else:
            not_exposed.append(candidate)
    return exposed, not_exposed


def filter_rhombus_by_exposure(planar, non_planar, atoms, height=2.5):
    """
    Given planar and non-planar candidate rhombus configurations, further filter them based on exposure.
    If a candidate is not exposed (i.e. is_exposed_rhombus() returns False), remove it from both groups
    and add it to the invalid_rhombus list.

    Parameters:
        planar (list): List of candidate rhombus configurations classified as planar.
        non_planar (list): List of candidate configurations classified as non-planar.
        atoms: An ASE Atoms object.
        height (float): Height threshold (in Å) for exposure.
    
    Returns:
        tuple: (planar_exposed, non_planar_exposed, invalid_rhombus)
            - planar_exposed: List of planar candidates that are exposed.
            - non_planar_exposed: List of non-planar candidates that are exposed.
            - invalid_rhombus: List of candidates that are not exposed.
    """
    exposed_planar, not_exposed_planar = filter_exposed_rhombus(planar, atoms, height)
    exposed_non_planar, not_exposed_non_planar = filter_exposed_rhombus(non_planar, atoms, height)
    invalid_rhombus = not_exposed_planar + not_exposed_non_planar
    return exposed_planar, exposed_non_planar, invalid_rhombus


def write_indices_to_file(file_path, indices_list):
    with open(file_path, 'w') as file:
        for indices in indices_list:
            text = ' '.join(map(str, indices))
            file.write(text + '\n')

def calculate_normal_vector_TOP(anchor_atom_index, cluster, cutoff):
    # Extract the position of the anchor atom
    anchor_atom = cluster[anchor_atom_index]
    anchor_position = anchor_atom.position
    
    # Find neighbors within the cutoff distance
    neighbors = []
    for neighbor in cluster:
        if not np.array_equal(neighbor.position, anchor_position):  # Ensure the neighbor is not the anchor atom
            vector = neighbor.position - anchor_position
            distance = np.linalg.norm(vector)
            if distance < cutoff:
                neighbors.append(neighbor.position)
    
    # Ensure we have enough neighbors to define a plane
    if len(neighbors) < 3:
        raise ValueError("Not enough neighbors to define a plane.")
    
    # Perform PCA on the neighbors to find the plane
    neighbors = np.array(neighbors)
    pca = PCA(n_components=3)
    pca.fit(neighbors)
    
    # The normal vector to the plane is the third principal component
    normal_vector = pca.components_[2]
    
    # Normalize the normal vector
    normal_vector /= np.linalg.norm(normal_vector)
    
    return normal_vector

def calculate_normal_vector_bri(atom_index, cluster, cutoff):
    vector_sum = np.zeros(3)
    atom = cluster[atom_index]

    for neighbor in cluster:
        if neighbor != atom:
            vector = atom.position - neighbor.position
            distance = np.linalg.norm(vector)
            if distance < cutoff and distance > 0:  # Avoid zero division
                vector_sum += vector / distance

    if np.linalg.norm(vector_sum) > 0:
        return vector_sum / np.linalg.norm(vector_sum)
    else:
        return np.array([0, 0, 1])  # Default to z-direction

def count_short_metal_nonmetal_bonds(cluster, metal='Ru', max_distance=1.6):
    """
    Count the number of bonds shorter than a specific distance between metal atoms and non-metal atoms in the cluster.

    Parameters:
    - cluster (ASE Atoms object): The cluster of atoms.
    - metal_symbol (str): Symbol of the metal atoms in the cluster.
    - max_distance (float): The maximum distance threshold for bond counting.

    Returns:
    - count (int): Number of bonds shorter than max_distance between metal atoms and non-metal atoms.
    """
    # Find indices of metal and non-metal atoms
    metal_indices = [i for i, atom in enumerate(cluster) if atom.symbol == metal]
    nonmetal_indices = [i for i, atom in enumerate(cluster) if atom.symbol != metal]

    count = 0

    # Calculate distances between metal and non-metal atoms
    for i in metal_indices:
        for j in nonmetal_indices:
            distance = cluster.get_distance(i, j, mic=True)
            if distance < max_distance:
                count += 1
    return count

def find_shortest_bond(atoms, cutoff=2.5):
    # read POSCAR
    # atoms = read(path)
    
    # get the index for N
    n_index = None
    ru_indices = []
    
    for i, atom in enumerate(atoms):
        if atom.symbol == 'N':
            n_index = i
        elif atom.symbol == 'Ru':
            ru_indices.append(i)

    if n_index is None:
        raise ValueError("No nitrogen (N) atom found in the POSCAR file.")  

    #Calculate the distance between N and each Ru atom, and filter out the distances that are less than the cutoff.
    distances = []

    for ru_index in ru_indices:
        distance = atoms.get_distance(n_index, ru_index, mic=True)
        if distance < cutoff:
            distances.append((ru_index, distance))

    
    ##"Sort by distance and find the shortest distance."
    if distances:
        distances.sort(key=lambda x: x[1])
        shortest_distance = distances[0]
        num_short_bonds = len(distances)
    else:
        shortest_distance = None
        num_short_bonds = 0

    return shortest_distance, num_short_bonds


##################=========================

def add_anchor_to_single_site_bri(poscar_path,exposed_indexs, metal='Ru', anchor_atom='N', output_prefix='POSCAR'):
    cluster = read(poscar_path)
    cutoff = 3.0  # Adjust based on your system
    modified_cluster = copy.deepcopy(cluster)
    def get_geo(modified_cluster, exposed_indexs):

        A = modified_cluster[exposed_indexs[0] ]  # Adjust index for 0-based indexing
        B = modified_cluster[exposed_indexs[1] ]  
        
        translation_vector = B.position - A.position
        # Calculate the normal vector for positioning the N atom
        normal_vector = calculate_normal_vector_bri(exposed_indexs[0] , modified_cluster, cutoff)
    
        # Position for N atom
        anchor_position = A.position + normal_vector * 1.6 +  0.5 * translation_vector  
        
        # Add N atom to the structure
        modified_cluster += Atoms(anchor_atom, positions=[anchor_position])
        return modified_cluster

    modified_cluster = copy.deepcopy(cluster)
    cluster_1 = get_geo(modified_cluster, exposed_indexs)
    num_bonds_1 = count_short_metal_nonmetal_bonds(cluster_1, metal, max_distance=1.6) 

    modified_cluster = copy.deepcopy(cluster)
    cluster_2 = get_geo(modified_cluster, exposed_indexs[::-1])
    
    num_bonds_2 = count_short_metal_nonmetal_bonds(cluster_2, metal, max_distance=1.6) 

    final_cluster = cluster_1
    if num_bonds_1 > num_bonds_2: 
        final_cluster = cluster_2 

    shortest_distance, num_short_bonds = find_shortest_bond(final_cluster, 2.5)
    # print(shortest_distance, num_short_bonds)
    if shortest_distance[1] > 1.5 and num_short_bonds < 3:     
        destination_path = os.path.join('bri', '_'.join([str(i) for i in exposed_indexs])) 
        os.makedirs(destination_path, exist_ok=True)
        output_file = destination_path + '/POSCAR'
    
        write(output_file, final_cluster, format='vasp', vasp5=True)

def add_one_bri(path, exposed_indices):
    poscar_path = path +'/POSCAR_bare'
    cluster = read(poscar_path)
    os.makedirs(path + '/bri', exist_ok=True)
    for combination in itertools.combinations(exposed_indices, 2):
        ## get the atoms that form the hollow site 
        sorted_combination = sorted(combination)
        bri_atom_indices = [i  for i in sorted_combination] # indice count from 0
        bri_atoms_positions = [cluster[i].position for i in bri_atom_indices]

        if not is_valid_bri(bri_atoms_positions):
            continue
        
        add_anchor_to_single_site_bri(poscar_path,sorted_combination,output_prefix='POSCAR')
 
def get_new_coordiates(cluster, A1, A2, B1, list_H):
    '''A1, A2, B1 are the coordinates, list_H is the list of coordinates '''
    # Calculate B2 as the center (midpoint) of C and D
    B2 = calculate_triangle_center(list_H) 
    translation_vector = A2 - B1
    B2_translated = B2 + translation_vector 
    
    # Step 2: Calculate vectors a and b
    a = A2 - A1
    a_normalized = a / np.linalg.norm(a)
    b = B2_translated - A2  # New b after moving B1 to A2
    b_normalized = b / np.linalg.norm(a)
    
    # Calculate rotation axis and angle
    axis = np.cross(b_normalized, a_normalized)
    
    axis_normalized = axis / np.linalg.norm(axis) if np.linalg.norm(axis) != 0 else b_normalized
    
    angle = np.arccos(np.clip(np.dot(b_normalized, a_normalized), -1.0, 1.0))     # Angle of rotation
    
    # Apply rotation around A2
    rotation = R.from_rotvec(axis_normalized * angle)
    
    new_list_H = [rotation.apply(i - B1) + A2 for i  in cluster.get_positions()[1:]]
    
    return new_list_H# C_final, D_final, E_final

def add_more_atoms(site):
    ''' Add species with more than three atoms to the anchor atoms. This function works for top, bri, and hollows sites '''
    sites = [int(i)  for i in site.split('_')]
    atoms = read(site + '/POSCAR')
    coordinates = [atoms[i].position for i in sites]
    A1 = calculate_triangle_center(coordinates)  # fcc site 
    A2 = atoms[-1].position     # N site, add N first and then put the NH,NH2,NH3 et al using the N site as reference.

    atoms_temp = read('../POSCAR_temp')
    B1 = atoms_temp[0].position

    # Initialize list_H to store atoms forming bonds with atoms_temp[0]
    list_H = []
    # Define bond distance threshold
    bond_distance = 1.9  # in Ångstroms
    positions = atoms_temp.get_positions()

    # Iterate over all atoms except the first one (atoms_temp[0])
    for i in range(1, len(positions)):
        atom_position = positions[i]
        distance = euclidean(B1, atom_position)
        
        if distance <= bond_distance:
            list_H.append(atom_position)
    
    # Convert list_H to numpy array if needed
    list_H = np.array(list_H)
    
    new_list_H  = get_new_coordiates(atoms_temp, A1, A2, B1, list_H)
    list_ele = atoms_temp.get_chemical_symbols()[1:]
    out_geo = copy.deepcopy(atoms)
    
    for num, coord in enumerate(new_list_H):
        out_geo += Atoms(list_ele[num], positions=[coord])
        
    write(site+'/POSCAR_bri', out_geo, format='vasp', vasp5=True)

def add_more_bri(path):
    """Add larger moelcules to all the bridge sites. The sites are the folders generated by the add_one_bri function"""
    os.chdir(path + '/bri')
    all_items = os.listdir('.')
    folders = [item for item in all_items if os.path.isdir(item)]
    folders = [folder for folder in folders if folder.count('_') == 1]
        
    for site in folders:     
        add_more_atoms(site)
        
def get_active_sites(cluster_path,metal = 'Ru'):   # Path: the cluster model directory
    poscar = os.path.join(cluster_path, 'POSCAR')
    atoms = read(poscar)
    bond_length_threshold = 2.7
    

    """ Step1: obtain the top sites""" 
    list_file = os.path.join(cluster_path, 'list')   ### list all the top sites in one line
    if not os.path.exists(list_file):
        # print("Warning: No list file found. Examine the exposed sites using Coordination Numbers.")
        top_indices = get_top_sites(atoms, metal, mult=0.9)  ## Seems that this method is not ideal
    else:
        with open(list_file, 'r') as f_in:
            top_indices = f_in.readline().rstrip().split()
            top_indices = [int(i) for i in top_indices]
            surface_indices = [i-1 for i in top_indices]
 
    """ Step2: obtain the bridge sites """     
    bridge_list = []
    for combination in itertools.combinations(top_indices, 2):
        ## get the atoms that form the hollow site 
        sorted_combination = sorted(combination)
        distance = atoms.get_distance(sorted_combination[0], sorted_combination[1])
        if distance < 2.6:
            bridge_list.append(sorted_combination)

    l_rhombus  = []
    for bridge_site in bridge_list: 
        B, D  =  bridge_site[:]
        indices = find_rhombus(atoms, top_indices, B, D, bond_threshold=3.0)
        l_rhombus.append(indices)
       
    l_rhombus = [item for item in l_rhombus if item] ## remove the empty 

    ############## Step3: categorize  rhombus sites: planar or edge
    planar, non_planar = filter_rhombus_by_dihedral(l_rhombus, atoms, dihedral_threshold=120)
    new_planar, new_non_planar, invalid_rhombus = filter_rhombus_by_AC_BD(planar, non_planar, atoms, AC_cutoff=3.0)
    # coplanar_threshold = 0.5  # Adjust threshold as needed, unit is Å
    # coplanar_rhombus = filter_coplanar_rhombuses(l_rhombus, atoms, coplanar_threshold)
    # edge_rhombus = filter_rhombuses_by_dihedral(l_rhombus, atoms, 40, 120)

    write_indices_to_file(os.path.join(cluster_path, 'all_sites.txt'), l_rhombus)
    write_indices_to_file(os.path.join(cluster_path, 'planar_sites.txt'), new_planar)
    write_indices_to_file(os.path.join(cluster_path, 'edge_sites.txt'), new_non_planar)
    write_indices_to_file(os.path.join(cluster_path, 'invalid_sites.txt'), invalid_rhombus)

    return l_rhombus, new_planar, new_non_planar



#### Construct adsorption configurations

# Define bonding criteria (in Angstrom)
R_N_BOND_MAX = 2.4  # Max Ru-N bond distance
R_H_BOND_MAX = 2.2  # Max Ru-H bond distance
N_N_BOND_MAX = 1.4  # Max N-N bond distance for N2 detection

def get_bonded_atoms(atoms, atom_index, element, cutoff):
    """Find all atoms of a given element bonded to a specified atom."""
    bonded = []
    for i, atom in enumerate(atoms):
        if atom.symbol == element and i != atom_index:
            dist = atoms.get_distance(atom_index, i)
            if dist < cutoff:
                bonded.append(i)
    return bonded



def classify_single_atom_adsorption(atoms, element, bond_cutoff):
    """Classify adsorption of single atoms (N or H) based on Ru coordination."""
    atom_indices = [i for i, atom in enumerate(atoms) if atom.symbol == element]
    results = []

    for idx in atom_indices:
        ru_bonded = get_bonded_atoms(atoms, idx, "Ru", bond_cutoff)

        # Determine adsorption type
        if len(ru_bonded) == 1:
            adsorption_type = "top"
        elif len(ru_bonded) == 2:
            adsorption_type = "bridge"
        elif len(ru_bonded) == 3:
            adsorption_type = "hollow"
        else:
            adsorption_type = "unknown"

        ru_site_str = format_ru_sites(ru_bonded)
        results.append((idx, adsorption_type, ru_site_str))

    # return results
    return ru_site_str


def format_ru_sites(ru_indices):
    """Convert a list of Ru indices into a sorted string format."""
    ru_indices = [i + 1 for i in ru_indices]  # Convert to 1-based index
    ru_indices.sort()  # Sort in ascending order
    return "_".join(map(str, ru_indices)) if ru_indices else "N/A"

def is_valid_bridge3(atoms, n1, n2, ru_bonded_n1, ru_bonded_n2):
    """
    Check if N2 forms a Bridge-3 adsorption:
    - Both N atoms bind to the same two Ru atoms.
    - N2 bond is at an angle between 75°-90° relative to Ru-Ru bond.

    Returns:
        True  -> Valid Bridge-3
        False -> Not a Bridge-3
    """
    if set(ru_bonded_n1) != set(ru_bonded_n2):
        return False  # Ensure both N atoms are bonded to the same two Ru atoms

    ru1, ru2 = ru_bonded_n1  # Get the two Ru atoms

    # Get atomic positions
    n1_pos = atoms[n1].position
    n2_pos = atoms[n2].position
    ru1_pos = atoms[ru1].position
    ru2_pos = atoms[ru2].position

    # Compute N2 and Ru-Ru vectors
    n2_vector = n2_pos - n1_pos
    n2_vector /= np.linalg.norm(n2_vector)

    ru_vector = ru2_pos - ru1_pos
    ru_vector /= np.linalg.norm(ru_vector)

    # Compute angle in degrees
    angle_rad = np.arccos(np.clip(np.dot(n2_vector, ru_vector), -1.0, 1.0))
    angle_deg = np.degrees(angle_rad)

    return 75 <= angle_deg <= 115  # Allow angle range between 75°-90°



def is_valid_fcc0(atoms, n1, n2, ru_bonded_n1, ru_bonded_n2):
    """
    Check if N2 forms an FCC-0 adsorption:
    - One N (N1 or N2) binds to two Ru atoms (A, B).
    - The other N (N2 or N1) binds to a third, separate Ru atom (C), instead of both A and B.

    Returns:
        True  -> Valid FCC-0
        False -> Not an FCC-0 site (likely Bridge-2)
    """
    if len(ru_bonded_n1) != 2 or len(ru_bonded_n2) != 1:
        if len(ru_bonded_n1) != 1 or len(ru_bonded_n2) != 2:
            return False  # Ensure N1 binds to two Ru, and N2 binds to one OR vice versa

    # Get the two Ru atoms bound to N1 (A, B) and the single Ru bound to N2 (C)
    ru_n1_set = set(ru_bonded_n1)
    ru_n2_set = set(ru_bonded_n2)

    # **For FCC-0: N2 must bind to a separate Ru (C) not in (A, B)**
    return len(ru_n1_set.intersection(ru_n2_set)) == 0

def is_valid_rhombus_site(atoms, n1, n2, ru_bonded_n1, ru_bonded_n2):
    """
    Check if N2 forms a Rhombus adsorption:
    - N1 binds to three Ru atoms (A, B, C).
    - N2 binds to one of (A, B, C) and a new Ru atom (D).
    - Distinguishes from FCC-2 where N2 binds to two of (A, B, C).

    Returns:
        True  -> Valid Rhombus site
        False -> Not a Rhombus site (likely FCC-2)
    """
    if len(ru_bonded_n1) != 3 or len(ru_bonded_n2) != 2:
        return False  # Ensure correct Ru coordination

    # Get the three Ru atoms bound to N1 (A, B, C)
    ru_n1_set = set(ru_bonded_n1)

    # Get the two Ru atoms bound to N2
    ru_n2_set = set(ru_bonded_n2)

    # **Check if N2 binds to only one of (A, B, C) and one new Ru (D)**
    common_rus = ru_n1_set.intersection(ru_n2_set)
    unique_rus = ru_n2_set.difference(ru_n1_set)

    return len(common_rus) == 1 and len(unique_rus) == 1



def classify_N2_adsorption(atoms):
    """Classify N2 adsorption types and format Ru sites as strings."""
    nitrogen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == "N"]
    # print('N2', nitrogen_indices)
    # Pair N atoms into N2 molecules based on N-N bonding distance
    n2_molecules = []
    visited = set()
    for i in nitrogen_indices:
        if i in visited:
            continue
        for j in nitrogen_indices:
            if i != j and atoms.get_distance(i, j) < N_N_BOND_MAX:
                n2_molecules.append((i, j))
                visited.add(i)
                visited.add(j)
                break

    results = []
    
    for n1, n2 in n2_molecules:
        
        ru_bonded_n1 = get_bonded_atoms(atoms, n1, "Ru", R_N_BOND_MAX)
        ru_bonded_n2 = get_bonded_atoms(atoms, n2, "Ru", R_N_BOND_MAX)
        
        # **Normalize Bond Counts** (Handle N1/N2 ordering issue)
        bond_counts = tuple(sorted([len(ru_bonded_n1), len(ru_bonded_n2)]))
        
        
        # **Check for Gas-Phase N2**
        if atoms.get_distance(n1, n2) < 1.3 and len(ru_bonded_n1) == 0 and len(ru_bonded_n2) == 0:
            adsorption_type = "gas"
            ru_site_str = "N/A"
        
        # **Classify Adsorption Type**
        elif bond_counts == (0, 1):
            adsorption_type = "top"
            ru_site_str = format_ru_sites(ru_bonded_n1 if len(ru_bonded_n1) == 1 else ru_bonded_n2)
        elif bond_counts == (1, 1) and ru_bonded_n1 == ru_bonded_n2:
            adsorption_type = "top"
            ru_site_str = format_ru_sites(ru_bonded_n1)  # One Ru index
        elif bond_counts == (1, 1):
            adsorption_type = "bridge-1"
            ru_site_str = format_ru_sites(list(set(ru_bonded_n1 + ru_bonded_n2)))  # Two Ru atoms
        
        elif bond_counts == (0, 2):
            # **Bridge-0: One N binds to two Ru atoms, the other is in the gas phase**
            adsorption_type = "bridge-0"
            ru_site_str = format_ru_sites(set(ru_bonded_n1 + ru_bonded_n2))
    
        elif bond_counts == (1, 2):
            # **Distinguish between FCC-0 and Bridge-2**
            if is_valid_fcc0(atoms, n1, n2, ru_bonded_n1, ru_bonded_n2):
                adsorption_type = "fcc-0"
            else:
                adsorption_type = "bridge-2"
            ru_site_str = format_ru_sites(set(ru_bonded_n1 + ru_bonded_n2))
        elif bond_counts == (2, 2):
            # **Bridge-3 detection (Both N bind to the same two Ru atoms & angle between 75°-90°)**
            if is_valid_bridge3(atoms, n1, n2, ru_bonded_n1, ru_bonded_n2):
                adsorption_type = "bridge-3"
                ru_site_str = format_ru_sites(set(ru_bonded_n1 + ru_bonded_n2))
            else:
                if len(set(ru_bonded_n1 + ru_bonded_n2)) == 3:
                    adsorption_type = "fcc-1"
                    ru_site_str = format_ru_sites(list(set(ru_bonded_n1 + ru_bonded_n2)))  # Three Ru atoms
                elif len(set(ru_bonded_n1 + ru_bonded_n2)) == 4:
                    adsorption_type = "trapezoid-1"
                    ru_site_str = format_ru_sites(list(set(ru_bonded_n1 + ru_bonded_n2)))  # Three Ru atoms
                    
                else:
                    adsorption_type = "unknown"
                    ru_site_str = "N/A"
        elif bond_counts == (2, 3):
            # **Distinguish between FCC-2 and Rhombus**
            if is_valid_rhombus_site(atoms, n1, n2, ru_bonded_n1, ru_bonded_n2):
                adsorption_type = "rhombus"
            else:
                adsorption_type = "fcc-2"
            ru_site_str = format_ru_sites(set(ru_bonded_n1 + ru_bonded_n2))  # Three Ru atoms
        elif bond_counts == (1, 3):
            adsorption_type = "fcc-3"
            ru_site_str = format_ru_sites(list(set(ru_bonded_n1 + ru_bonded_n2)))  # Three Ru atoms
        else:
            adsorption_type = "unknown"
            ru_site_str = "N/A"

        #results.append((n1, n2, adsorption_type, ru_site_str))
        results.append((adsorption_type, ru_site_str))

    return results   


def get_shortest_ru_n_distance(atoms):
    """Get the shortest Ru-N distance in the given structure."""
    ru_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'Ru']
    n_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'N']

    min_distance = float('inf')
    for ru in ru_indices:
        for n in n_indices:
            distance = atoms.get_distance(ru, n)
            if distance < min_distance:
                min_distance = distance

    return min_distance

def add_N2_top_sites(atoms, n_ru_distance=2.0, n_n_distance=1.19, metal='Ru'):
    """
    Add N2 molecules at top adsorption sites (Top1 and Top2).
    """

    connections, cn_of_connected_atoms, top_sites, bridge_sites, hollow_sites = get_connection(atoms, metal=metal, mult=0.9)
    mass_center = np.mean([atom.position for atom in atoms if atom.symbol == 'Ru'], axis=0)

    for i, ru_index in enumerate(top_sites):
        pos_top = atoms[ru_index].position  # Ru atom serving as the top site
        direction = pos_top - mass_center
        direction /= np.linalg.norm(direction)

        if np.dot(direction, pos_top - mass_center) < 0:
            direction = -direction  # Flip direction if needed

        ### **Top1 Placement**
        n1_position_top1 = pos_top + direction * n_ru_distance
        n2_position_top1 = n1_position_top1 + direction * n_n_distance
        n2_top1 = Atoms('N2', positions=[n1_position_top1, n2_position_top1])

        new_atoms_top1 = atoms.copy()
        new_atoms_top1.extend(n2_top1)

        if get_shortest_ru_n_distance(new_atoms_top1) < 1.5:
            write(f"POSCAR_top1_{i+1}_check1", new_atoms_top1, format='vasp')
            direction = -direction
            n1_position_top1 = pos_top + direction * n_ru_distance
            n2_position_top1 = n1_position_top1 + direction * n_n_distance
            n2_top1 = Atoms('N2', positions=[n1_position_top1, n2_position_top1])
            new_atoms_top1 = atoms.copy()
            new_atoms_top1.extend(n2_top1)

            if get_shortest_ru_n_distance(new_atoms_top1) < 1.5:
                write(f"POSCAR_top1_{i+1}_check2", new_atoms_top1, format='vasp')
            else:
                write(f"POSCAR_top1_{i+1}.vasp", new_atoms_top1, format='vasp')
        else:
            write(f"POSCAR_top1_{i+1}.vasp", new_atoms_top1, format='vasp')

        ### **Top2 Placement**
        nn_center = pos_top + direction * n_ru_distance
        perp_vector = np.cross(direction, [0, 0, 1])
        if np.linalg.norm(perp_vector) < 1e-6:
            perp_vector = np.cross(direction, [1, 0, 0])
        perp_vector /= np.linalg.norm(perp_vector)

        n1_position_top2 = nn_center - perp_vector * (n_n_distance / 2)
        n2_position_top2 = nn_center + perp_vector * (n_n_distance / 2)
        n2_top2 = Atoms('N2', positions=[n1_position_top2, n2_position_top2])

        new_atoms_top2 = atoms.copy()
        new_atoms_top2.extend(n2_top2)

        write(f"POSCAR_top2_{i+1}.vasp", new_atoms_top2, format='vasp')

def add_N2_bridge_sites(atoms, n_ru_distance=2.15, n_n_distance=1.19):
    """
    Add N2 molecules at bridge adsorption sites (Bridge-1, Bridge-2, Bridge-3, FCC-0).
    """

    connections, cn_of_connected_atoms, top_sites, bridge_sites, hollow_sites = get_connection(atoms, metal=metal, mult=0.9)
    mass_center = np.mean([atom.position for atom in atoms if atom.symbol == 'Ru'], axis=0)

    for i, (ru1, ru2) in enumerate(bridge_sites):
        pos1, pos2 = atoms[ru1].position, atoms[ru2].position
        bridge_mid = (pos1 + pos2) / 2
        direction = bridge_mid - mass_center
        direction /= np.linalg.norm(direction)

        ### **Bridge-1 Placement (Parallel to Ru-Ru)**
        nn_center = bridge_mid + direction * n_ru_distance
        ru_vector = pos2 - pos1
        ru_vector /= np.linalg.norm(ru_vector)

        n1_position_bridge1 = nn_center - ru_vector * (n_n_distance / 2)
        n2_position_bridge1 = nn_center + ru_vector * (n_n_distance / 2)
        n2_bridge1 = Atoms('N2', positions=[n1_position_bridge1, n2_position_bridge1])

        new_atoms_bridge1 = atoms.copy()
        new_atoms_bridge1.extend(n2_bridge1)

        if get_shortest_ru_n_distance(new_atoms_bridge1) < 1.5:
            write(f"POSCAR_bridge1_{i+1}_check1", new_atoms_bridge1, format='vasp')
            direction = -direction
            n1_position_bridge1 = nn_center - ru_vector * (n_n_distance / 2)
            n2_position_bridge1 = nn_center + ru_vector * (n_n_distance / 2)
            n2_bridge1 = Atoms('N2', positions=[n1_position_bridge1, n2_position_bridge1])
            new_atoms_bridge1 = atoms.copy()
            new_atoms_bridge1.extend(n2_bridge1)

            if get_shortest_ru_n_distance(new_atoms_bridge1) < 1.5:
                write(f"POSCAR_bridge1_{i+1}_check2", new_atoms_bridge1, format='vasp')
            else:
                write(f"POSCAR_bridge1_{i+1}.vasp", new_atoms_bridge1, format='vasp')
        else:
            write(f"POSCAR_bridge1_{i+1}.vasp", new_atoms_bridge1, format='vasp')

        ### **Bridge-2 & FCC-0 Placement**
        possible_c_rus = [ru for ru in connections[ru1] if ru != ru2] + [ru for ru in connections[ru2] if ru != ru1]

        if possible_c_rus:  # Ensure a third Ru atom exists
            ru3 = possible_c_rus[0]  # Select one possible C atom
            pos3 = atoms[ru3].position

            if is_valid_fcc0(atoms, ru1, ru2, [ru1, ru2], [ru3]):  
                adsorption_type = "fcc-0"
                filename = f"POSCAR_fcc0_{i+1}.vasp"
            else:
                adsorption_type = "bridge-2"
                filename = f"POSCAR_bridge2_{i+1}.vasp"

            n1_position = (pos1 + pos2) / 2 + direction * n_ru_distance
            n2_position = pos3 + direction * n_ru_distance  # N2 binds to Ru3 in FCC-0 or Ru1/Ru2 in Bridge-2

            n2_new = Atoms('N2', positions=[n1_position, n2_position])

            new_atoms = atoms.copy()
            new_atoms.extend(n2_new)

            if get_shortest_ru_n_distance(new_atoms) < 1.5:
                write(f"{filename.replace('.vasp', '_check1.vasp')}", new_atoms, format='vasp')
                direction = -direction
                n1_position = (pos1 + pos2) / 2 + direction * n_ru_distance
                n2_position = pos3 + direction * n_ru_distance
                n2_new = Atoms('N2', positions=[n1_position, n2_position])

                new_atoms = atoms.copy()
                new_atoms.extend(n2_new)

                if get_shortest_ru_n_distance(new_atoms) < 1.5:
                    write(f"{filename.replace('.vasp', '_check2.vasp')}", new_atoms, format='vasp')
                else:
                    write(filename, new_atoms, format='vasp')
            else:
                write(filename, new_atoms, format='vasp')

        ### **Bridge-3 Placement (N₂ perpendicular to Ru-Ru, both N have 2 Ru-N bonds)**
        perp_vector = np.cross(ru_vector, direction)
        perp_vector /= np.linalg.norm(perp_vector)

        nn_center_bridge3 = bridge_mid + direction * n_ru_distance
        n1_position_bridge3 = nn_center_bridge3 + perp_vector * (n_n_distance / 2)
        n2_position_bridge3 = nn_center_bridge3 - perp_vector * (n_n_distance / 2)
        n2_bridge3 = Atoms('N2', positions=[n1_position_bridge3, n2_position_bridge3])

        new_atoms_bridge3 = atoms.copy()
        new_atoms_bridge3.extend(n2_bridge3)

        if get_shortest_ru_n_distance(new_atoms_bridge3) < 1.5:
            write(f"POSCAR_bridge3_{i+1}_check1", new_atoms_bridge3, format='vasp')
            direction = -direction
            n1_position_bridge3 = nn_center_bridge3 + perp_vector * (n_n_distance / 2)
            n2_position_bridge3 = nn_center_bridge3 - perp_vector * (n_n_distance / 2)
            n2_bridge3 = Atoms('N2', positions=[n1_position_bridge3, n2_position_bridge3])
            new_atoms_bridge3 = atoms.copy()
            new_atoms_bridge3.extend(n2_bridge3)

            if get_shortest_ru_n_distance(new_atoms_bridge3) < 1.5:
                write(f"POSCAR_bridge3_{i+1}_check2", new_atoms_bridge3, format='vasp')
            else:
                write(f"POSCAR_bridge3_{i+1}.vasp", new_atoms_bridge3, format='vasp')
        else:
            write(f"POSCAR_bridge3_{i+1}.vasp", new_atoms_bridge3, format='vasp')

        print(f"Saved POSCAR_bridge3_{i+1}")
  

def add_hollow_sites(atoms, n_ru_distance=1.95, n_n_distance=1.19, n1_height=1.5,):
    """
    Add N2 molecules at hollow adsorption sites (FCC-1, FCC-2, FCC-3).
    """

    connections, cn_of_connected_atoms, top_sites, bridge_sites, hollow_sites = get_connection(atoms, metal=metal, mult=0.9)
    ru_positions = np.array([atom.position for atom in atoms if atom.symbol == 'Ru'])
    mass_center = np.mean(ru_positions, axis=0)

    for i, (ru1, ru2, ru3) in enumerate(hollow_sites):
        triangle_sites = [ru1, ru2, ru3]

        if not is_exposed_triangle(atoms, triangle_sites, height=2.5):
            print(f"Skipping hollow site {ru1+1}, {ru2+1}, {ru3+1} (not exposed)")
            continue

        pos1, pos2, pos3 = atoms[ru1].position, atoms[ru2].position, atoms[ru3].position
        hollow_center = np.mean([pos1, pos2, pos3], axis=0)

        # **Compute outward direction from mass center to hollow site**
        direction = hollow_center - mass_center
        direction /= np.linalg.norm(direction)

        if np.dot(direction, hollow_center - mass_center) < 0:
            direction = -direction  # Flip direction vector if needed

        normal = np.cross(pos2 - pos1, pos3 - pos1)
        normal /= np.linalg.norm(normal)

        ### **FCC-1 Placement**
        n1_position_fcc1 = (pos1 + pos2) / 2 + direction * n_ru_distance
        n2_position_fcc1 = (pos2 + pos3) / 2 + direction * n_ru_distance
        n2_fcc1 = Atoms('N2', positions=[n1_position_fcc1, n2_position_fcc1])

        new_atoms_fcc1 = atoms.copy()
        new_atoms_fcc1.extend(n2_fcc1)

        if get_shortest_ru_n_distance(new_atoms_fcc1) < 1.5:
            write(f"POSCAR_fcc1_{i+1}_check1", new_atoms_fcc1, format='vasp')
            direction = -direction
            n1_position_fcc1 = (pos1 + pos2) / 2 + direction * n_ru_distance
            n2_position_fcc1 = (pos2 + pos3) / 2 + direction * n_ru_distance
            n2_fcc1 = Atoms('N2', positions=[n1_position_fcc1, n2_position_fcc1])
            new_atoms_fcc1 = atoms.copy()
            new_atoms_fcc1.extend(n2_fcc1)

            if get_shortest_ru_n_distance(new_atoms_fcc1) < 1.5:
                write(f"POSCAR_fcc1_{i+1}_check2", new_atoms_fcc1, format='vasp')
            else:
                write(f"POSCAR_fcc1_{i+1}.vasp", new_atoms_fcc1, format='vasp')
        else:
            write(f"POSCAR_fcc1_{i+1}.vasp", new_atoms_fcc1, format='vasp')

        ### **FCC-2 Placement**
        n1_position_fcc2 = hollow_center + normal * n1_height
        perp_vector = np.cross(normal, direction)
        perp_vector /= np.linalg.norm(perp_vector)
        n2_position_fcc2 = n1_position_fcc2 + perp_vector * n_n_distance
        n2_fcc2 = Atoms('N2', positions=[n1_position_fcc2, n2_position_fcc2])

        new_atoms_fcc2 = atoms.copy()
        new_atoms_fcc2.extend(n2_fcc2)

        if get_shortest_ru_n_distance(new_atoms_fcc2) < 1.5:
            write(f"POSCAR_fcc2_{i+1}_check1", new_atoms_fcc2, format='vasp')
            direction = -direction
            n2_position_fcc2 = n1_position_fcc2 + perp_vector * n_n_distance
            n2_fcc2 = Atoms('N2', positions=[n1_position_fcc2, n2_position_fcc2])
            new_atoms_fcc2 = atoms.copy()
            new_atoms_fcc2.extend(n2_fcc2)

            if get_shortest_ru_n_distance(new_atoms_fcc2) < 1.5:
                write(f"POSCAR_fcc2_{i+1}_check2", new_atoms_fcc2, format='vasp')
            else:
                write(f"POSCAR_fcc2_{i+1}.vasp", new_atoms_fcc2, format='vasp')
        else:
            write(f"POSCAR_fcc2_{i+1}.vasp", new_atoms_fcc2, format='vasp')

        ### **FCC-3 Placement**
        n1_position_fcc3 = hollow_center + normal * n1_height
        closest_ru = min([pos1, pos2, pos3], key=lambda p: np.linalg.norm(p - n1_position_fcc3))
        n2_position_fcc3 = n1_position_fcc3 + perp_vector * n_n_distance
        n2_fcc3 = Atoms('N2', positions=[n1_position_fcc3, n2_position_fcc3])

        new_atoms_fcc3 = atoms.copy()
        new_atoms_fcc3.extend(n2_fcc3)

        if get_shortest_ru_n_distance(new_atoms_fcc3) < 1.5:
            write(f"POSCAR_fcc3_{i+1}_check1", new_atoms_fcc3, format='vasp')
            direction = -direction
            n2_position_fcc3 = n1_position_fcc3 + perp_vector * n_n_distance
            n2_fcc3 = Atoms('N2', positions=[n1_position_fcc3, n2_position_fcc3])
            new_atoms_fcc3 = atoms.copy()
            new_atoms_fcc3.extend(n2_fcc3)

            if get_shortest_ru_n_distance(new_atoms_fcc3) < 1.5:
                write(f"POSCAR_fcc3_{i+1}_check2", new_atoms_fcc3, format='vasp')
            else:
                write(f"POSCAR_fcc3_{i+1}.vasp", new_atoms_fcc3, format='vasp')
        else:
            write(f"POSCAR_fcc3_{i+1}.vasp", new_atoms_fcc3, format='vasp')
