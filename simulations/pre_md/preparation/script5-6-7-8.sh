
gmx=/usr/local/gromacs-5.1.2_plumed-2.3b/bin/gmx_mpi


mkdir 5-equil
cd 5-equil

cp ../param_files/equil.mdp .
cp ../4-min.ions/start.top .
cp ../4-min.ions/confout.gro .
cp ../4-min.ions/*.itp .
$gmx genrestr -f confout.gro -o posre.itp <<eof
1
eof

$gmx grompp -f equil.mdp -c confout.gro -p start.top -o equil.tpr -maxwarn 2

$gmx mdrun -s equil.tpr -v -ntomp 6

rm ener.edr
rm traj.trr
rm mdout.mdp
rm traj.xtc
rm state.cpt
rm state_prev.cpt

cd ..

mkdir 6-nvt 
cd 6-nvt

cp ../param_files/pre.nvt.mdp .
cp ../5-equil/start.top .
cp ../5-equil/confout.gro .
cp ../5-equil/*.itp .

$gmx grompp -f pre.nvt.mdp -c confout.gro -p start.top -o nvt.tpr -maxwarn 2
$gmx mdrun -s nvt.tpr -v -ntomp 6

rm ener.edr
rm traj.trr
rm mdout.mdp
rm traj.xtc
rm state.cpt
rm state_prev.cpt

cd ..

mkdir 7-npt
cd 7-npt

cp ../param_files/pre.npt.mdp .
cp ../6-nvt/start.top .
cp ../6-nvt/confout.gro .
cp ../6-nvt/*.itp .

$gmx grompp -f pre.npt.mdp -c confout.gro -p start.top -o npt.tpr -maxwarn 2
$gmx mdrun -s npt.tpr -v -ntomp 6

rm ener.edr
rm traj.trr
rm mdout.mdp
rm traj.xtc
rm state.cpt
rm state_prev.cpt


cd ..

mkdir 8-eq
cd 8-eq

cp ../param_files/eq.mdp .
cp ../7-npt/start.top .
cp ../7-npt/confout.gro .
cp ../7-npt/*.itp .

$gmx grompp -f eq.mdp -c confout.gro -p start.top -o eq.tpr -maxwarn 2

$gmx mdrun -s eq.tpr -v -ntomp 6 

$gmx trjconv -f confout.gro -s eq.tpr -ur compact -pbc mol -o confout.pdb <<eof
0
eof

#rm ener.edr
#rm md.log
#rm traj.trr
#rm mdout.mdp
#rm traj.xtc
#rm state.cpt
#rm state_prev.cpt



