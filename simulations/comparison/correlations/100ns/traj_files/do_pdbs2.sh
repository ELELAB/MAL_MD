labels=( 9-md_1000ns_2 )

for l in ${labels[@]}; do 
	echo "20" |  gmx_mpi trjconv -s ${l}.tpr -f ${l}.xtc -o ${l}_A.frame0.pdb -dump 0 -n index.ndx 
        echo "22" |  gmx_mpi trjconv -s ${l}.tpr -f ${l}.xtc -o ${l}_B.frame0.pdb -dump 0 -n index.ndx
done
