#!/bin/bash

filename=ions
seed=101
##########################


mkdir 3-ions
cd 3-ions

while [ $seed -lt 400 ];
do
cp ../param_files/ions.mdp
cp ../param_files/script3-ions.py
cp ../2-solvate/start.gro .
cp ../2-solvate/start.top .
cp ../2-solvate/*.itp .
grompp -f ions.mdp -c start.gro -p start.top -o ions.tpr

genion -s ions.tpr  -neutral -conc 0.1 -seed $seed -o subions_$seed.gro -p start.top  << eof
15
eof

grompp -f ions.mdp -c subions_$seed.gro -p start.top -o ions_ion.tpr

trjconv -s ions_ion.tpr -f subions_$seed.gro -pbc mol -ur compact -o subions_$seed.pdb << eof
0
eof

pymol -c subions_$seed.pdb script3-ions.py

if [ $? -eq 0 ];
then
echo "really bad________________________________________________________________________________"
rm subions_$seed.pdb
rm subions_$seed.gro
rm genion.log
rm start.gro
rm *.top ions_ion.tpr
rm \#*
seed=$(( $seed + 1 ))
else
seed=$(( $seed + 2000 ))
rm \#*
echo "I finally did it!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
pymol subions_*.pdb script.pml
fi
done
