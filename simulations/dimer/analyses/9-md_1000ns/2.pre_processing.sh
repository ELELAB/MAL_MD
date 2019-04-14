#!/bin/bash

path_gmx=/usr/local/gromacs-5.1.2_plumed-2.3b/bin/gmx_mpi

xtc=traj_centered
tpr=sim
gro=confout

$path_gmx make_ndx -f $gro.gro -o prot_MG.ndx << eof
1 | 13
keep 20
q
eof

$path_gmx make_ndx -f $gro.gro -o prot_MG_A.ndx << eof
ri 1-414
keep 20
q
eof

$path_gmx make_ndx -f $gro.gro -o prot_MG_B.ndx << eof
ri 415-828
keep 20
q
eof



$path_gmx trjconv -f $xtc.xtc -s $tpr.tpr -pbc mol -ur compact -o traj_centered_prot_MG_A.xtc -n prot_MG_A.ndx
$path_gmx trjconv -f $xtc.xtc -s $tpr.tpr -pbc mol -ur compact -o traj_centered_prot_MG_B.xtc -n prot_MG_B.ndx 

$path_gmx convert-tpr -s $tpr.tpr -n prot_MG_A.ndx -o sim_prot_MG_A.tpr
$path_gmx convert-tpr -s $tpr.tpr -n prot_MG_B.ndx -o sim_prot_MG_B.tpr


$path_gmx trjconv -f traj_centered_prot_MG_A.xtc -s sim_prot_MG_A.tpr -fit rot+trans -dt 1500 -o movie.dt1500_A.xtc <<eof
3
1
eof
$path_gmx trjconv -f traj_centered_prot_MG_B.xtc -s sim_prot_MG_B.tpr -fit rot+trans -dt 1500 -o movie.dt1500_B.xtc <<eof
3
1
eof

