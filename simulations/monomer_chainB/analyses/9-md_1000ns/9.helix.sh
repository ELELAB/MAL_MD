#!/bin/bash

path_gmx=/usr/local/gromacs-5.1.2_plumed-2.3b/bin/gmx_mpi

xtc=traj_centered_prot_MG
tpr=sim_prot_MG

mkdir 9.helix
cd 9.helix

$path_gmx make_ndx -f ../2.*/$tpr.tpr -o<<eof
q
eof

$path_gmx helix -f ../2.*/$xtc.xtc -s ../2.*/$tpr.tpr -n index.ndx<<eof
1
eof 



