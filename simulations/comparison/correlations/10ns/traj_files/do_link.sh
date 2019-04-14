systems=(9-md_1000ns  9-md_1000ns_2 9-md_1000ns_3)
labels=(9-md_1000ns 9-md_1000ns_2 9-md_1000ns_3)
xtc=traj_centered.xtc
tpr=sim.tpr

for s in $(seq 0 3); do 
        unlink  ${labels[$s]}.xtc
        unlink  ${labels[$s]}.tpr
	ln -s ../../../../dimer/analyses/${systems[$s]}/$xtc ${labels[$s]}.xtc
	ln -s ../../../../dimer/analyses/${systems[$s]}/$tpr ${labels[$s]}.tpr
done
