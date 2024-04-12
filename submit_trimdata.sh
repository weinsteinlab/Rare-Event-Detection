#!/bin/bash
#SBATCH -p edison,panda_physbio
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G

systems="ligand"
start_traj=replica
further_strdie=1

scriptfolder='/home/hex4001/RED'
python3 $scriptfolder/trimdata.py $systems $start_traj $further_strdie
exit
