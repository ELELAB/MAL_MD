#!/bin/bash

path_gmx=/usr/local/gromacs-5.1.2_plumed-2.3b/bin/gmx_mpi

mkdir dihedrals
cd dihedrals
$path_gmx chi -f ../2.*/traj_centered_prot_MG.xtc -s ../2.*/sim_prot_MG.tpr  -phi -psi -omega -all -maxchi 6

mkdir omega psi phi chi
mv omega*.xvg omega
mv histo-omega* omega
mv chi*.xvg chi
mv histo-chi*.xvg chi
mv psi*.xvg psi
mv histo-psi*.xvg psi
mv phi*.xvg phi
mv histo-phi*.xvg phi

#python chi_histo_plot.py -f config.cfg
