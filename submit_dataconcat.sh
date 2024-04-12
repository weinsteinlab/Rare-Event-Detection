#!/bin/bash
#SBATCH -p edison,panda_physbio
#SBATCH --job-name=input3_input4_
#SBATCH --output=input3_input4.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G

scriptfolder='/home/hex4001/RED'
python3 $scriptfolder/dataconcat.py input1 input2 input3 input4
#python3 dataconcat_chol.py input1 input2 input3 input4
#python3 dataconcat_chol2.py input1 input2 input3 input4
exit
