#!/bin/bash
labels=( 9-md_1000ns_2 )

for l in ${labels[@]}; do 
	echo '20' | gmx_mpi trjconv -s ${l}.tpr -f ${l}.xtc -o ${l}_A.xtc -pbc whole -dt 100 -n index.ndx
	echo '21' | gmx_mpi trjconv -s ${l}.tpr -f ${l}.xtc -o ${l}_A.frame0.CA.pdb -dump 0 -n index.ndx
        python add_segname -n A -p ${l}_A.frame0.CA.pdb
        echo '22' | gmx_mpi trjconv -s ${l}.tpr -f ${l}.xtc -o ${l}_B.xtc -pbc whole -dt 100 -n index.ndx
        echo '23' | gmx_mpi trjconv -s ${l}.tpr -f ${l}.xtc -o ${l}_B.frame0.CA.pdb -dump 0 -n index.ndx 
	python add_segname -n B -p ${l}_B.frame0.CA.pdb
done
