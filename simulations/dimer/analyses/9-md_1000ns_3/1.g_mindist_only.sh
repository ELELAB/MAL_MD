#!/bin/bash

path_gmx=/usr/local/gromacs-5.1.2_plumed-2.3b/bin/gmx_mpi
#path_files=/data/user/shared_projects/p53_mutants/pre_MD_2xwr_S215G_NO_DNA/md_analysis/9-md
path_analysis=/data/user/shared_projects/chalmers/mal/simulations/dimer/analyses/9-md_1000ns_2/
path_index=/data/user/shared_projects/chalmers/mal/simulations/dimer/analyses/9-md_1000ns_2/
path_plot=/data/user/shared_projects/p53_mutants/PLOT/standard-plot-PCA/

j=traj_comp
k=sim


$path_gmx mindist -f traj_centered.xtc -s $k.tpr -od mindist_center.xvg -pi << eof
1
eof

