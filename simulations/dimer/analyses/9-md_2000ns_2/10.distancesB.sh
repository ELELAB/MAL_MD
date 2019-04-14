mkdir dist_B
cd dist_B
gmx_mpi make_ndx -f ../sim_prot_MG_B.tpr -o index_distB.ndx <<eof
ri 745
ri 487  
ri 494
q
eof


gmx_mpi pairdist -f ../traj_centered_prot_MG_B.xtc -s ../sim_prot_MG_B.tpr -n index_distB.ndx -ref -sel -refgrouping res -selgrouping res -o dist_K331.Q73_K331.R80.xvg
