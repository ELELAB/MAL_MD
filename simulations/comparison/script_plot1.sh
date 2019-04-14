cd PCA_domains



cd PCA_Prot_1 
sed '/^#/ d' proj2d_1vs2.xvg > parz
sed '/^@/ d' parz > proj2d_1vs2_m.xvg
cp /data/user/shared_projects/chalmers/mal/simulations/comparison/ref_time.xvg .
paste ref_time.xvg proj2d_1vs2_m.xvg | awk '{print $1,$3,$4}' > proj2d_1vs2_mt.xvg
rm parz
sed '/^#/ d' proj2d_1vs3.xvg > parz
sed '/^@/ d' parz > proj2d_1vs3_m.xvg
paste ref_time.xvg proj2d_1vs3_m.xvg | awk '{print $1,$3,$4}' > proj2d_1vs3_mt.xvg
rm parz *.eps

cp ../../plot_pca1v2_single_comp.pnl .
cp ../../plot_pca1v3_single_comp.pnl .

gnuplot plot_pca1v2_single_comp.pnl
gnuplot plot_pca1v3_single_comp.pnl


cd  ../PCA_cterm_dom_6 
sed '/^#/ d' proj2d_1vs2.xvg > parz
sed '/^@/ d' parz > proj2d_1vs2_m.xvg
cp /data/user/shared_projects/chalmers/mal/simulations/comparison_pca/ref_time.xvg .
paste ref_time.xvg proj2d_1vs2_m.xvg | awk '{print $1,$3,$4}' > proj2d_1vs2_mt.xvg
rm parz
sed '/^#/ d' proj2d_1vs3.xvg > parz
sed '/^@/ d' parz > proj2d_1vs3_m.xvg
paste ref_time.xvg proj2d_1vs3_m.xvg | awk '{print $1,$3,$4}' > proj2d_1vs3_mt.xvg
rm parz *.eps

cp ../../plot_pca1v2_single_comp.pnl .
cp ../../plot_pca1v3_single_comp.pnl .

gnuplot plot_pca1v2_single_comp.pnl
gnuplot plot_pca1v3_single_comp.pnl

 
cd  ../PCA_loop_long_2 
sed '/^#/ d' proj2d_1vs2.xvg > parz
sed '/^@/ d' parz > proj2d_1vs2_m.xvg
cp /data/user/shared_projects/chalmers/mal/simulations/comparison_pca/ref_time.xvg .
paste ref_time.xvg proj2d_1vs2_m.xvg | awk '{print $1,$3,$4}' > proj2d_1vs2_mt.xvg
rm parz
sed '/^#/ d' proj2d_1vs3.xvg > parz
sed '/^@/ d' parz > proj2d_1vs3_m.xvg
paste ref_time.xvg proj2d_1vs3_m.xvg | awk '{print $1,$3,$4}' > proj2d_1vs3_mt.xvg
rm parz *.eps

cp ../../plot_pca1v2_single_comp.pnl .
cp ../../plot_pca1v3_single_comp.pnl .

gnuplot plot_pca1v2_single_comp.pnl
gnuplot plot_pca1v3_single_comp.pnl


cd  ../PCA_loop_long_shortdef_3  
sed '/^#/ d' proj2d_1vs2.xvg > parz
sed '/^@/ d' parz > proj2d_1vs2_m.xvg
cp /data/user/shared_projects/chalmers/mal/simulations/comparison_pca/ref_time.xvg .
paste ref_time.xvg proj2d_1vs2_m.xvg | awk '{print $1,$3,$4}' > proj2d_1vs2_mt.xvg
rm parz
sed '/^#/ d' proj2d_1vs3.xvg > parz
sed '/^@/ d' parz > proj2d_1vs3_m.xvg
paste ref_time.xvg proj2d_1vs3_m.xvg | awk '{print $1,$3,$4}' > proj2d_1vs3_mt.xvg
rm parz *.eps

cp ../../plot_pca1v2_single_comp.pnl .
cp ../../plot_pca1v3_single_comp.pnl .

gnuplot plot_pca1v2_single_comp.pnl
gnuplot plot_pca1v3_single_comp.pnl

cd  ../PCA_loop_short_4  
sed '/^#/ d' proj2d_1vs2.xvg > parz
sed '/^@/ d' parz > proj2d_1vs2_m.xvg
cp /data/user/shared_projects/chalmers/mal/simulations/comparison_pca/ref_time.xvg .
paste ref_time.xvg proj2d_1vs2_m.xvg | awk '{print $1,$3,$4}' > proj2d_1vs2_mt.xvg
rm parz
sed '/^#/ d' proj2d_1vs3.xvg > parz
sed '/^@/ d' parz > proj2d_1vs3_m.xvg
paste ref_time.xvg proj2d_1vs3_m.xvg | awk '{print $1,$3,$4}' > proj2d_1vs3_mt.xvg
rm parz *.eps

cp ../../plot_pca1v2_single_comp.pnl .
cp ../../plot_pca1v3_single_comp.pnl .

gnuplot plot_pca1v2_single_comp.pnl
gnuplot plot_pca1v3_single_comp.pnl


cd  ../PCA_nterm_dom_5  
sed '/^#/ d' proj2d_1vs2.xvg > parz
sed '/^@/ d' parz > proj2d_1vs2_m.xvg
cp /data/user/shared_projects/chalmers/mal/simulations/comparison_pca/ref_time.xvg .
paste ref_time.xvg proj2d_1vs2_m.xvg | awk '{print $1,$3,$4}' > proj2d_1vs2_mt.xvg
rm parz
sed '/^#/ d' proj2d_1vs3.xvg > parz
sed '/^@/ d' parz > proj2d_1vs3_m.xvg
paste ref_time.xvg proj2d_1vs3_m.xvg | awk '{print $1,$3,$4}' > proj2d_1vs3_mt.xvg
rm parz *.eps

cp ../../plot_pca1v2_single_comp.pnl .
cp ../../plot_pca1v3_single_comp.pnl .

gnuplot plot_pca1v2_single_comp.pnl
gnuplot plot_pca1v3_single_comp.pnl


cd  ../PCA_Prot_Core_1
sed '/^#/ d' proj2d_1vs2.xvg > parz
sed '/^@/ d' parz > proj2d_1vs2_m.xvg
cp /data/user/shared_projects/chalmers/mal/simulations/comparison_pca/ref_time.xvg .
paste ref_time.xvg proj2d_1vs2_m.xvg | awk '{print $1,$3,$4}' > proj2d_1vs2_mt.xvg
rm parz
sed '/^#/ d' proj2d_1vs3.xvg > parz
sed '/^@/ d' parz > proj2d_1vs3_m.xvg
paste ref_time.xvg proj2d_1vs3_m.xvg | awk '{print $1,$3,$4}' > proj2d_1vs3_mt.xvg
rm parz *.eps

cp ../../plot_pca1v2_single_comp.pnl .
cp ../../plot_pca1v3_single_comp.pnl .

gnuplot plot_pca1v2_single_comp.pnl
gnuplot plot_pca1v3_single_comp.pnl
