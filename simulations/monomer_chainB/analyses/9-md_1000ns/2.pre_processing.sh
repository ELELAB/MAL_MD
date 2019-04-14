#!/bin/bash

path_gmx=/usr/local/gromacs-5.1.2_plumed-2.3b/bin/gmx_mpi

xtc=traj_centered
tpr=sim
gro=confout

mkdir pre_processing
cd pre_processing


$path_gmx make_ndx -f ../$gro.gro -o prot_MG.ndx << eof
1 | 13
keep 20
q
eof


$path_gmx trjconv -f ../1.*/$xtc.xtc -s ../$tpr.tpr -pbc mol -ur compact -o traj_centered_prot_MG.xtc -n prot_MG.ndx 

$path_gmx convert-tpr -s ../$tpr.tpr -n prot_MG.ndx -o sim_prot_MG.tpr


$path_gmx trjconv -f traj_centered_prot_MG.xtc -s sim_prot_MG.tpr -fit rot+trans -dt 2000 -o movie_prot_MG.dt2000.pdb -n prot_MG.ndx

#mkdir models
#cd models
#trjconv -f ../traj_mol_ur_nojump_PROT.xtc -s sim_PROT.tpr -fit rot+trans -dt 10000 -o model.pdb -sep -n chains.ndx << eof
#28
#eof



