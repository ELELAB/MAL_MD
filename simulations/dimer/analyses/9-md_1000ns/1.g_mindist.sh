#!/bin/bash

path_gmx=/usr/local/gromacs-5.1.2_plumed-2.3b/bin/gmx_mpi
path_plot=/data/user/shared_projects/p53_mutants/PLOT/standard-plot-PCA/

j=traj_comp
k=sim

$path_gmx make_ndx -f $k.tpr -o center_atom.ndx << eof
r 356 & 3
q
eof

$path_gmx trjconv -center -f $j.xtc -o traj_centered.xtc -n center_atom.ndx -s $k.tpr -pbc mol -ur tric << eof
20
0
eof

$path_gmx mindist -f traj_centered.xtc -s $k.tpr -od mindist_center.xvg -pi << eof
1
eof


