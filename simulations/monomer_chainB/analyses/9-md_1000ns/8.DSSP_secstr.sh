#!/bin/bash
export DSSP=/usr/local/bin/dssp-2.2.8

path_gmx=/usr/local/gromacs-5.1.2_plumed-2.3b/bin/gmx_mpi

xtc=traj_centered_prot_MG
tpr=sim_prot_MG
output=dssp_monomerB

mkdir dssp
cd dssp


cp ../dssp_eps_parameters.m2p .

#File parameters.m2p for .eps was taken from: (and modified)
#http://ringo.ams.sunysb.edu/index.php/MD_Simulation:_Protein_in_Water_(Pt._2)

$path_gmx do_dssp -f ../2.*/$xtc.xtc -s ../2.*/$tpr.tpr -sc scount_dt100.xvg -dt 100 -o $output_dt100.xpm  <<eof
1
eof 

$path_gmx xpm2ps -f $output.xpm -di dssp_eps_parameters.m2p -o $output.eps

convert eps:$output.eps png:$output.png


