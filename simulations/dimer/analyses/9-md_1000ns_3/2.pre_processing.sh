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
#trjconv -f $xtc.xtc -s $tpr.tpr -pbc mol -ur compact -o traj_mol_ur.xtc -n prot_Mg.ndx

#tpbconv -s $tpr.tpr -n prot_Mg.ndx -o sim_prot_Mg.tpr

#trjconv -f traj_mol_ur.xtc -s sim_prot_Mg.tpr -pbc nojump -o traj_mol_ur_nojump.xtc -n prot_Mg.ndx



$path_gmx trjconv -f $xtc.xtc -s $tpr.tpr -pbc mol -ur compact -o traj_centered_prot_MG_A.xtc -n prot_MG_A.ndx
$path_gmx trjconv -f $xtc.xtc -s $tpr.tpr -pbc mol -ur compact -o traj_centered_prot_MG_B.xtc -n prot_MG_B.ndx 

$path_gmx convert-tpr -s $tpr.tpr -n prot_MG_A.ndx -o sim_prot_MG_A.tpr
$path_gmx convert-tpr -s $tpr.tpr -n prot_MG_B.ndx -o sim_prot_MG_B.tpr

#trjconv -f traj_mol_ur_PROT.xtc -s sim_PROT.tpr -pbc nojump -o traj_mol_ur_nojump_PROT.xtc -n prot.ndx 

#rm -f traj_mol_ur.xtc traj_mol_ur_PROT.xtc

#g_rms -f traj_mol_ur_nojump.xtc -s sim_prot_Mg.tpr -o rmsd.xvg << eof
#5
#5
#eof
#g_rms -f traj_mol_ur_nojump.xtc -s sim_prot_Mg.tpr -n chains.ndx -o rmsd_chainA.xvg <<eof
#3
#26
#eof
#g_rms -f traj_mol_ur_nojump.xtc -s sim_prot_Mg.tpr -n chains.ndx -o rmsd_chainB.xvg <<eof
#3
#27
#eof
#g_gyrate -f traj_mol_ur_nojump.xtc -s sim_prot_Mg.tpr -o gyrate.xvg <<eof
#1
#eof


$path_gmx trjconv -f traj_centered_prot_MG_A.xtc -s sim_prot_MG_A.tpr -fit rot+trans -dt 1500 -o movie.dt1500_A.xtc <<eof
3
1
eof
$path_gmx trjconv -f traj_centered_prot_MG_B.xtc -s sim_prot_MG_B.tpr -fit rot+trans -dt 1500 -o movie.dt1500_B.xtc <<eof
3
1
eof

#mkdir models
#cd models
#trjconv -f ../traj_mol_ur_nojump_PROT.xtc -s sim_PROT.tpr -fit rot+trans -dt 10000 -o model.pdb -sep -n chains.ndx << eof
#28
#eof



