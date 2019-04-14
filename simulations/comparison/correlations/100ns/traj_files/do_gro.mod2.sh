labels=( 9-md_1000ns_2 )

for l in ${labels[@]}; do 
        echo "0" | gmx_mpi trjconv -s ${l}.tpr -f ${l}.xtc -o ${l}.frame0.gro -n index.ndx -dump 0
done
