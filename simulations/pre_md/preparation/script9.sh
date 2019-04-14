
gmx=/usr/local/gromacs-5.1.2_plumed-2.3b/bin/gmx_mpi

mkdir 9-md
cd 9-md

cp ../param_files/sim.mdp .
cp ../8-eq/start.top .
cp ../8-eq/confout.gro .
cp ../8-eq/*.itp .

$gmx grompp -f sim.mdp -c confout.gro -p start.top -o sim.tpr -maxwarn 2

cd ..

