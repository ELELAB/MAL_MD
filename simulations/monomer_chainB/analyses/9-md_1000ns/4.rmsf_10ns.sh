#!/bin/bash

path_gmx=/usr/local/gromacs-5.1.2_plumed-2.3b/bin/gmx_mpi


mkdir rmsf
cd rmsf
j=0
z=10000
for ((i=1; i <=100; i++)); do
$path_gmx rmsf -f ../2.pre_processing/traj_centered_prot_MG.xtc -s ../2.pre_processing/sim_prot_MG.tpr -res -o rmsf$i.xvg -b $j -e $z << eof
1
eof
sed -i '/#/d' rmsf$i.xvg
sed -i '/@/d' rmsf$i.xvg
mv rmsf$i.xvg 00$i.rmsf.oct
j=$(($j+10000))
z=$(($z+10000))
done
awk 'FNR==1{f++}{a[f,FNR]=$2}END{for(x=1;x<=FNR;x++){for(y=1;y<ARGC;y++)printf("%s ",a[y,x]);print ""}}' 0*.rmsf.oct > sum.rmsf.oct
paste 001.rmsf.oct sum.rmsf.oct > all.rmsf.oct
cat all.rmsf.oct  | awk '{ s = 0; for (i = 2; i <= NF; i++) s += $i; print $1, (NF > 1) ? s / (NF - 1) : 0; }' > AV_RMSF_chainB.oct
#paste rmsf1.oct rmsf2.oct rmsf3.oct rmsf4.oct rmsf5.oct rmsf6.oct rmsf7.oct rmsf8.oct rmsf9.oct rmsf10.oct | awk '{print $1, " ",$2, " ",$4, " ",$6, " ", $8, " ", $10, " ", $12, " ", $14, " ", $16, " ", $18, " ", $20}' > all.rmsf.oct
#rm all.rmsf.oct
