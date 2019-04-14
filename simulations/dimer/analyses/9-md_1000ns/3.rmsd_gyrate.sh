#!/bin/bash

path_gmx=/usr/local/gromacs-5.1.2_plumed-2.3b/bin/gmx_mpi

xtc=traj_centered
tpr=sim
gro=confout


mkdir rmsd_gyrate
cd rmsd_gyrate

$path_gmx make_ndx -f ../sim_prot_MG_A.tpr -o index_core_A.ndx <<eof
r 177-187 | r 209-225 | r 235-238 | r 251-265 | r 270-273 | r 281-297 | r 303-306 | r 313-322 | r 327-330 | r 338-351 | r 355-357 | r 365-378 | r 382-384
14 & 5
q 
eof

$path_gmx make_ndx -f ../sim_prot_MG_B.tpr -o index_core_B.ndx <<eof
r 177-187 | r 209-225 | r 235-238 | r 251-265 | r 270-273 | r 281-297 | r 303-306 | r 313-322 | r 327-330 | r 338-351 | r 355-357 | r 365-378 | r 382-384
14 & 5
q 
eof


$path_gmx rms -f ../traj_centered_prot_MG_A.xtc -s ../sim_prot_MG_A.tpr -o rmsd_mainchain_A.xvg -n index_core_A.ndx << eof
5
5
eof
$path_gmx rms -f ../traj_centered_prot_MG_A.xtc -s ../sim_prot_MG_A.tpr -o rmsd_core_A.xvg -n index_core_A.ndx << eof
15
5
eof


$path_gmx rms -f ../traj_centered_prot_MG_A.xtc -s ../sim_prot_MG_A.tpr -o rmsd_MG_A.xvg -n index_core_A.ndx <<eof
15
13
eof

$path_gmx gyrate -f ../traj_centered_prot_MG_A.xtc -s ../sim_prot_MG_A.tpr -o gyrate_A.xvg <<eof
1
eof


$path_gmx rms -f ../traj_centered_prot_MG_B.xtc -s ../sim_prot_MG_B.tpr -o rmsd_mainchain_B.xvg -n index_core_B.ndx << eof
5
5
eof
$path_gmx rms -f ../traj_centered_prot_MG_B.xtc -s ../sim_prot_MG_B.tpr -o rmsd_core_B.xvg -n index_core_B.ndx << eof
15
5
eof


$path_gmx rms -f ../traj_centered_prot_MG_B.xtc -s ../sim_prot_MG_B.tpr -o rmsd_MG_B.xvg -n index_core_B.ndx <<eof
15
13
eof

$path_gmx gyrate -f ../traj_centered_prot_MG_B.xtc -s ../sim_prot_MG_B.tpr -o gyrate_B.xvg <<eof
1
eof

cd ..
mkdir models
cd models
$path_gmx trjconv -f ../traj_centered_prot_MG_A.xtc -s ../sim_prot_MG_A.tpr -fit rot+trans -dt 10000 -o model_A.pdb -sep -n ../prot_MG_A.ndx
$path_gmx trjconv -f ../traj_centered_prot_MG_B.xtc -s ../sim_prot_MG_B.tpr -fit rot+trans -dt 10000 -o model_B.pdb -sep -n ../prot_MG_B.ndx 



