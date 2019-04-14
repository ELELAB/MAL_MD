#!/bin/bash

path_gmx=/usr/local/gromacs-5.1.2_plumed-2.3b/bin/gmx_mpi

xtc=traj_centered
tpr=sim
gro=confout

mkdir rmsd_gyrate
cd rmsd_gyrate

#$path_gmx trjconv -f $xtc.xtc -s $tpr.tpr -pbc mol -ur compact -o traj_centered_prot_MG.xtc -n prot_MG.ndx 

#$path_gmx convert-tpr -s $tpr.tpr -n prot_MG.ndx -o sim_prot_MG.tpr

$path_gmx make_ndx -f ../$tpr.tpr -o index_core.ndx <<eof
ri 177-187 | ri 209-225 | ri 235-238 | ri 251-265 | ri 270-273 | ri 281-297 | ri 303-306 | ri 313-322 | ri 327-330 | ri 338-351 | ri 355-357 | ri 365-378 | ri 382-384
20 & 5
q 
eof

$path_gmx rms -f ../2.*/traj_centered_prot_MG.xtc -s ../2.*/sim_prot_MG.tpr -o rmsd_mainchain.xvg -n index_core.ndx << eof
5
5
eof
$path_gmx rms -f ../2.*/traj_centered_prot_MG.xtc -s ../2.*/sim_prot_MG.tpr -o rmsd_core.xvg -n index_core.ndx << eof
21
5
eof


$path_gmx rms -f ../2.*/traj_centered_prot_MG.xtc -s ../2.*/sim_prot_MG.tpr -o rmsd_MG.xvg -n index_core.ndx <<eof
21
13
eof

$path_gmx gyrate -f ../2.*/traj_centered_prot_MG.xtc -s ../2.*/sim_prot_MG.tpr -o gyrate.xvg <<eof
1
eof

mkdir models
cd models

$path_gmx trjconv -f ../../2.*/traj_centered_prot_MG.xtc -s ../../2.*/sim_prot_MG.tpr -fit rot+trans -dt 10000 -o model.pdb -sep -n ../../2.*/prot_MG.ndx 



