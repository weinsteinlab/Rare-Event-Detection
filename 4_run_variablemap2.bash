#!/bin/bash
#SBATCH -p edison,panda_physbio
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G

systems="CHLp2b2"
start_traj=0
end_traj=11

scriptfolder='/home/hex4001/RED'
python3 $scriptfolder/variablemap2.py $systems $start_traj $end_traj
exit
