#!/bin/bash
#SBATCH -p edison,panda_physbio
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G

systems="ligand"
start_traj=replica
scriptfolder='/home/hex4001/RED'

python3 $scriptfolder/variablemap1.py $systems $start_traj
exit
