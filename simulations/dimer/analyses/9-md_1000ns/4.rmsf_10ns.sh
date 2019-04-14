#!/bin/bash

path_gmx=/usr/local/gromacs-5.1.2_plumed-2.3b/bin/gmx_mpi


mkdir rmsf_A rmsf_B
cd rmsf_A
j=0
z=10000
for ((i=1; i <=100; i++)); do
$path_gmx rmsf -f ../traj_centered_prot_MG_A.xtc -s ../sim_prot_MG_A.tpr -res -o rmsf$i.xvg -b $j -e $z << eof
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
cat all.rmsf.oct  | awk '{ s = 0; for (i = 2; i <= NF; i++) s += $i; print $1, (NF > 1) ? s / (NF - 1) : 0; }' > AV_RMSF_chainA.oct
rm 00*.rmsf.oct


cd ../rmsf_B
j=0
z=10000
for ((i=1; i <=100; i++)); do
$path_gmx rmsf -f ../traj_centered_prot_MG_B.xtc -s ../sim_prot_MG_B.tpr -res -o rmsf$i.xvg -b $j -e $z << eof
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
rm 00*.rmsf.oct
