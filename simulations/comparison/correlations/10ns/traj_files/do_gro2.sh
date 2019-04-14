labels=( 9-md_1000ns_2 )

for l in ${labels[@]}; do 
	echo "20" | gmx_mpi trjconv -s ${l}.tpr -f ${l}.xtc -o ${l}_A.frame0.gro -n index.ndx -dump 0 
        echo "22" | gmx_mpi trjconv -s ${l}.tpr -f ${l}.xtc -o ${l}_B.frame0.gro -n index.ndx -dump 0
done
