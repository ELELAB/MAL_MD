export DSSP=/usr/local/bin/dssp-2.2.8

mkdir dssp_A
cd dssp_A
gro=confout


gmx_mpi do_dssp -s ../sim_prot_MG_A.tpr -f ../traj_centered_prot_MG_A.xtc -map -sc scount_dt100.xvg -dt 100 -o dssp_analsys_dt100 <<eof
1
eof

cd ..


mkdir dssp_B
cd dssp_B

gmx_mpi do_dssp -s ../sim_prot_MG_B.tpr -f ../traj_centered_prot_MG_B.xtc -map -sc scount_dt100.xvg -dt 100 -o dssp_analsys_dt100 <<eof
1
eof

cd ..


