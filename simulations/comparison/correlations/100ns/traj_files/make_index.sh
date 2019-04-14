gmx_mpi make_ndx -f ../../../../dimer/analyses/9-md_1000ns/sim.tpr  -o index.ndx <<eof
ri 1-413
20 & 3
ri 415-827
22 & 3
ri 209-384
24 & 3
ri 623-798
26 & 3
ri 209-384 | ri 623-798
28 & 3
q
eof

# 20 r_1-413             :  6364 atoms
# 21 r_1-413_&_C-alpha   :   413 atoms
# 22 r_415-827           :  6364 atoms
# 23 r_415-827_&_C-alpha :   413 atoms
# 24 r_209-384           :  2745 atoms
# 25 r_209-384_&_C-alpha :   176 atoms
# 26 r_623-798           :  2745 atoms
# 27 r_623-798_&_C-alpha :   176 atoms
# 28 r_209-384_r_623-798 :  5490 atoms
# 29 r_209-384_r_623-798_&_C-alpha:   352 atoms
