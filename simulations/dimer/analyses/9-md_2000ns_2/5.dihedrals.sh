#!/bin/bash

path_gmx=/usr/local/gromacs-5.1.2_plumed-2.3b/bin/gmx_mpi

mkdir dihedrals_A dihedrals_B
cd dihedrals_A
$path_gmx chi -f ../traj_centered_prot_MG_A.xtc -s ../sim_prot_MG_A.tpr  -phi -psi -omega -all -maxchi 6
cp ../config.cfg .
mkdir omega psi phi chi
mv omega*.xvg omega
mv chi*.xvg chi
mv psi*.xvg psi
mv phi*.xvg phi
python chi_histo_plot.py -f config.cfg

cd ../dihedrals_B
$path_gmx chi -f ../traj_centered_prot_MG_B.xtc -s ../sim_prot_MG_B.tpr  -phi -psi -omega -all -maxchi 6
cp ../config.cfg .
mkdir omega psi phi chi
mv omega*.xvg omega
mv chi*.xvg chi
mv psi*.xvg psi
mv phi*.xvg phi
python chi_histo_plot.py -f config.cfg


