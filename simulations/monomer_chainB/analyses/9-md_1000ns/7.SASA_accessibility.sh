#!/bin/bash

path_gmx=/usr/local/gromacs-5.1.2_plumed-2.3b/bin/gmx_mpi

xtc=traj_centered_prot_MG
tpr=sim_prot_MG

mkdir sasa
cd sasa


$path_gmx make_ndx -f ../2.*/$tpr.tpr -o index_catalitic.ndx <<eof
ri 73
ri 170
ri 172
ri 194
ri 329
ri 331
ri 356
ri 360
ri 361
ri 384
ri 73 | ri 170 | ri 172 | ri 194 | ri 329 | ri 331 | ri 356 | ri 360 | ri 361 | ri 384
q
eof

$path_gmx sasa -f ../2.*/traj_centered_prot_MG.xtc -s ../2.*/sim_prot_MG.tpr -n index_catalitic.ndx -o sasa.xvg -or sasa_res.xvg -surface -output <<eof
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



