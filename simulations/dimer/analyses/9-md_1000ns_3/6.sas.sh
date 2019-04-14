mkdir sas_A sas_B
cd sas_A
gmx_mpi make_ndx -f ../sim_prot_MG_A.tpr -o index_cat_A.ndx <<eof
r 73
r 170
r 172
r 194
r 329
r 331
r 356
r 360
r 361
r 384
14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23
q
eof

gmx_mpi sasa -f ../traj_centered_prot_MG_A.xtc -n index_cat_A.ndx -s ../sim_prot_MG_A.tpr -o sasa.xvg -or sasa_res.xvg -surface -output <<eof
1
14
15
16
17
18
19
20
21
22
23
24
1
eof

cd ..
cd sas_B
gmx_mpi make_ndx -f ../sim_prot_MG_B.tpr -o index_cat_B.ndx <<eof
r 73
r 170
r 172
r 194
r 329
r 331
r 356
r 360
r 361
r 384
14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23
q
eof

gmx_mpi sasa -f ../traj_centered_prot_MG_B.xtc -n index_cat_B.ndx -s ../sim_prot_MG_B.tpr -o sasa.xvg -or sasa_res.xvg -surface -output <<eof
1
14
15
16
17
18
19
20
21
22
23
24
1
eof




