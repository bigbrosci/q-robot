for i in CHG* WAVE* AE* e.* e_*  o.*  o_* err.*  out.* REPORT* PROCAR* PCDAT* p4vasp.log; do 
    for j in $(find . -name $i); do 
        rm $j 
    done 
done 
