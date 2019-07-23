#!/bin/bash

path_gmx=/usr/local/gromacs-5.1.2_plumed-2.3b/bin/gmx_mpi
#path_files=/data/user/shared_projects/p53_mutants/pre_MD_2xwr_S215G_NO_DNA/md_analysis/9-md
path_analysis=/data/user/shared_projects/chalmers/mal/simulations/comparison/
path_index=/data/user/shared_projects/chalmers/mal/simulations/comparison/
path_plot=/data/user/shared_projects/chalmers/mal/simulations/comparison/

tpr=sim_prot_MG_A.tpr
traj=traj_monA_monB_dim1A_dim1B_dim2A_dim2B_dim3A_dim3B_1kko1A_1kko1B.xtc


#cd $path_analysis
#$path_gmx make_ndx -f $path_index/$tpr -o index_domains.ndx <<eof
#del 2-100
#r 177-187 | r 209-225 | r 235-238 | r 251-265 | r 270-273 | r 281-297 | r 303-306 | r 313-322 | r 327-330 | r 338-351 | r 355-357 | r 365-378 | r 382-384
#r 12-51
#r 18-47
#r 70-85
#r 4-146
#r 150-410
#q 
#eof

mkdir PCA_domains_last_long_1kko
cd PCA_domains_last_long_1kko

mkdir PCA_Prot_1 PCA_Prot_Core_2 PCA_loop_long_3 PCA_loop_long_shortdef_4 PCA_loop_short_5 PCA_nterm_dom_6 PCA_cterm_dom_7


J=1
for ((i=1; i <=7; i++)); do

cd *_$i
 
$path_gmx covar -f $path_index/$traj -s $path_index/$tpr  -o eigenval.xvg -v eigenvec.trr -av average.pdb -l covar.log -ascii covar.dat -xpm covar.xpm -xpma covara.xpmi -mwa -n $path_analysis/index_domains.ndx << eof
2
$J
eof

$path_gmx anaeig -f $path_index/$traj -s $path_index/$tpr -v eigenvec.trr -2d proj2d_1vs2.xvg -first 1 -last 2 -n $path_analysis/index_domains.ndx << eof
2
$J
eof
$path_gmx anaeig -f $path_index/$traj -s $path_index/$tpr -v eigenvec.trr -2d proj2d_1vs3.xvg -first 1 -last 3 -n $path_analysis/index_domains.ndx << eof
2
$J
eof

sed '/^#/ d' proj2d_1vs2.xvg > parz
sed '/^@/ d' parz > proj2d_1vs2_m.xvg
cp /data/user/shared_projects/chalmers/mal/simulations/comparison/ref_time.xvg .
paste ref_time.xvg proj2d_1vs2_m.xvg | awk '{print $1,$6,$7}' > proj2d_1vs2_mt.xvg
rm parz
sed '/^#/ d' proj2d_1vs3.xvg > parz
sed '/^@/ d' parz > proj2d_1vs3_m.xvg
paste ref_time.xvg proj2d_1vs3_m.xvg | awk '{print $1,$6,$7}' > proj2d_1vs3_mt.xvg
rm parz

cp ../../plot_pca1v2_single_comp.pnl .
cp ../../plot_pca1v3_single_comp.pnl .
cp ../../plot_proj2d_1vs2.pnl .
cp ../../plot_proj2d_1vs3.pnl .


gnuplot plot_pca1v2_single_comp.pnl
gnuplot plot_pca1v3_single_comp.pnl



gnuplot plot_proj2d_1vs2.pnl
gnuplot plot_proj2d_1vs3.pnl

mkdir ED_eig1 ED_eig2 ED_eig3
cd ED_eig1
cp cp ../../../dividi_pdb.pl .
$path_gmx anaeig -f $path_index/$traj -s $path_index/$tpr  -v ../eigenvec.trr  -extr models.pdb  -first 1 -last 1 -nframes 10 -n  $path_analysis/index_domains.ndx << eof
2
$J
eof
perl dividi_pdb.pl
cd ..
cd ED_eig2
cp cp ../../../dividi_pdb.pl .
$path_gmx anaeig -f  $path_index/$traj -s $path_index/$tpr  -v ../eigenvec.trr  -extr models.pdb  -first 2 -last 2 -nframes 10 -n  $path_analysis/index_domains.ndx << eof
2
$J
eof
perl dividi_pdb.pl
cd ..

cd ED_eig3
cp cp ../../../dividi_pdb.pl .
$path_gmx anaeig -f $path_index/$traj -s $path_index/$tpr  -v ../eigenvec.trr  -extr models.pdb  -first 3 -last 3 -nframes 10 -n  $path_analysis/index_domains.ndx << eof
2
$J
eof
perl dividi_pdb.pl

cd ..
cd ..

J=$(($J+1))
done

