gmx=/usr/local/gromacs-5.1.2_plumed-2.3b/bin/gmx_mpi

mkdir 4-min.ions
cd 4-min.ions

cp ../min.mdp .
cp ../3-ions/*.itp .
cp ../3-ions/start.top .
cp ../3-ions/confout.gro .
cp ../3-ions/subions_* .


$gmx grompp -f min.mdp -c subions_*.gro -p start.top -o min.tpr -maxwarn 2

$gmx mdrun -s min.tpr -v -ntomp 6


rm ener.edr
rm traj.trr
rm mdout.mdp
rm #confout.gro.1#

cd ..



