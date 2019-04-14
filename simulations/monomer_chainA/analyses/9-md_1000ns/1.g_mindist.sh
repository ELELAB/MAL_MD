#!/bin/bash

path_gmx=/usr/local/gromacs-5.1.2_plumed-2.3b/bin/gmx_mpi
#path_files=/data/user/shared_projects/p53_mutants/pre_MD_2xwr_S215G_NO_DNA/md_analysis/9-md
path_analysis=/data/user/shared_projects/chalmers/mal/simulations/monomer_chainA/nosubs/analyses/9-md_1micro_OK/pca/
path_index=/data/user/shared_projects/chalmers/mal/simulations/monomer_chainA/nosubs/analyses/9-md_1micro_OK/
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

"""
xvg2octave min-dist-pbc_center.xvg
awk '$2 <= 0.9' min-dist-pbc_center.oct > minbad_pbc.oct
awk 'END { print NR }' minbad_pbc.oct > minbad_number_pbc.oct
rm minbad_pbc.oct
"""

