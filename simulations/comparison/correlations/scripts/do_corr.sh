w
!/bin/bash
source ~/.GM467G

#  File: do_corr.sh
#
#  Copyright (C) 2012 Matteo Tiberti <matteo.tiberti@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.

# We now need to calculate the DCCM/LMI using wordom.
# As is usually done with g_corr or m_covar, the correlation matrix will be calculated as an average matrix. We will have to slice the complete trajectory in shorter sections and calculate a single DCCM for each one. The following variable define the length in ps of the time windows.

partsl=10000
# OPTIONAL: filter xtc and regenerate the pdb to keep only the atoms you need. This greatly speeds up g_covar. Set filtergroup accordingly.
#filtergroup=18
#echo $filtergroup | trjconv -n -s topol_prot_eq5.tpr -f traj_prot_centered_fitted.xtc -o traj_prot_centered_fitted_filtered.xtc
#echo $filtergroup | trjconv -n -s topol_prot_eq5.tpr -f traj_prot_centered_1.pdb -o traj_prot_centered_filtered_1.pdb
#./add_segname -p traj_prot_centered_filtered_1.pdb -r reference.pdb
# as before you may use "-n AA" for single-chain systems.
# 6: generate single trjs. You may want to discard the last part of the trajectory, so that you have exactly m windows of length $partsl. To do this you can add the "-e" option to the next line. Otherwise, you may need to calculated a weighted mean using the trajectory lengths as weights (this is not supported at the moment).

#labels=(1NJQ_c22starzn 1NJQ_c27zn 2L1O_c22starzn 2L1O_c27zn)
#labels=(1NJQ_c27)
labels=(1NJQ_c27Li 2L1O_c27Li)


rm \#*

for l in ${labels[@]}; do 

rm -rf ${l}
mkdir ${l}
cd ${l}

#rm *.corr *.dat traj_cf* *.ana average*pdb traj_fitted.xtc

ln -s ../add_segname . 
ln -s ../corr2dat .
ln -s ../corr_avg .

wordom=wordom_0.22-rc3.x86_64
pdb=../trajs_tops/${l}_frame0.CACH3.pdb
xtc=../trajs_tops/${l}.CACH3.xtc
ndx=../index.ndx

echo -e "10\n0" | trjconv -s $pdb -f $xtc -fit rot+trans -o traj_fitted.xtc -n $ndx
echo "0" | trjconv -s $pdb -f traj_fitted.xtc -o traj_cf_part.xtc -split $partsl
# 7: for each one:
trjsn=$(ls traj_cf_part*.xtc|wc -l)
for i in $(seq 0 $(($trjsn-2)) ); do
        xtcpart=traj_cf_part$i.xtc
	xtcpartfit=traj_cf_part_fit$i.xtc
        dcdpart=traj_cf_part$i.dcd
# 7a: create a wordom input file
# notice the selection (--SELE), which involves C, CA and N atoms of each residue. As LEVEL is RES and MASS is 0, the non mass-weighted average position between these atoms will be considered. To use CAs only, you may use:
# --SELE /*/*/CA
        cat <<EOF > CORR_part$i.ana
BEGIN corr
--TITLE corr_part$i.corr
--TYPE LMI
--SELE /*/*/@(CA|CH3)
--MASS 0
END
EOF
        # 7b: generate average structures
        #echo -e "0" | g_covar -s $pdb -f $xtcpart -av average"$i".pdb -nofit; rm eigen* covar.log
	echo -e "10\n0" | trjconv -f $xtcpart -s $pdb -fit rot+trans -o $xtcpartfit -n $ndx
	#echo "0" | g_covar -s $pdb -f $xtcpartfit -av average"$i".pdb -nofit
	echo "0" | g_rmsf -s $pdb -f $xtcpartfit -ox average"$i".pdb -nofit; rm rmsf.xvg;
        ./add_segname -p average"$i".pdb -o average"$i".pdb -n A
        # 7c: convert xtc to dcd
        $wordom -imol average"$i".pdb -itrj $xtcpartfit -otrj $dcdpart -conv -nopbc
        # 7d calculate correlations
        $wordom -iA CORR_part"$i".ana -imol average"$i".pdb -itrj $dcdpart -nopbc
done
# 8: calculate average matrix
./corr_avg corr_part*.corr > corr-average.corr
# If you wish you can convert this matrix to the .dat format using the corr2dat script.
# Done! Proceed to phase 3
./corr2dat -c corr-average.corr -o corr-average.dat

cd ..
done
