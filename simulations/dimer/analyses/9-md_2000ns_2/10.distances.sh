

mkdir dist_A 
cd dist_A
gmx_mpi make_ndx -f ../sim_prot_MG_A.tpr -o index_distA.ndx <<eof
ri 331
ri 73  
ri 80
q
eof


gmx_mpi pairdist -f ../traj_centered_prot_MG_A.xtc -s ../sim_prot_MG_A.tpr -n index_distA.ndx -ref -sel -refgrouping res -selgrouping res -o dist_K331.Q73_K331.R80.xvg <<eof
14
15
16
q
eof

