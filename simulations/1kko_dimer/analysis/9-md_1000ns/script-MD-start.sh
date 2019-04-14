#!/bin/bash
#SBATCH -J MAL_UNBOUND 
#SBATCH -A 2018-3-331
#SBATCH -N 1
#SBATCH --ntasks-per-node=32    # 32 tasks on 1 node
#SBATCH -t 24:00:00        	# time limits: 24 hours
#SBATCH --error myJob.err      # std-error file
#SBATCH --output myJob.out     # std-output file
export OMP_NUM_THREADS=1
module swap PrgEnv-cray PrgEnv-gnu
module add gromacs/5.1.2

t_now=$(date +%d%m-%H%M)".1node-4px1ntomp.out"

aprun -n 32 gmx_mpi mdrun  -s sim.tpr -v -cpi -append -ntomp 1 -maxh 23.99 >& $t_now
