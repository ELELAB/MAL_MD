labels=(1NJQ_c22starzn 1NJQ_c27zn 1NJQ_c27Li 2L1O_c22starzn 2L1O_c27zn 2L1O_c27Li)

for l in ${labels[@]}; do 
	echo '2' | gmx trjconv -s ${l}.tpr -f ${l}.xtc -o ${l}.noH.xtc -pbc whole -dt 100
	echo '2' | gmx trjconv -s ${l}.tpr -f ${l}.xtc -o ${l}_frame0.noH.pdb -dump 0 
done
