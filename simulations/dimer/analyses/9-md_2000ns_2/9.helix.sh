#!/bin/bash

path_gmx=/usr/local/gromacs-5.1.2_plumed-2.3b/bin/gmx_mpi


mkdir helix_A
cd helix_A

$path_gmx make_ndx -f ../sim_prot_MG_A.tpr -o<<eof
q
eof

$path_gmx helix -f ../traj_centered_prot_MG_A.xtc -s ../sim_prot_MG_A.tpr -n index.ndx<<eof
1
eof 

cd ..

mkdir helix_B
cd helix_B

$path_gmx make_ndx -f ../sim.tpr -o<<eof
ri 415-828
keep 20
eof

$path_gmx helix -f ../traj_centered.xtc -s ../sim.tpr -n index.ndx<<eof
eof 



