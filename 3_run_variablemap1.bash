#!/bin/bash

## declare an array variable
folders='
CHLp2b2
'
scriptfolder='/home/hex4001/RED'
## now loop through the above array
for i in $folders; do
   for replica in `seq 0 11`; do
   replica0=`printf %04d $replica`
      echo "$i $replica"
      cp $scriptfolder/submit_variablemap1.sh temp_${i}_${replica0}.sh
      sed -i "s/ligand/$i/g" temp_${i}_${replica0}.sh
      sed -i "s/replica/$replica0/g" temp_${i}_${replica0}.sh
      sbatch temp_${i}_${replica0}.sh
      rm temp_${i}_${replica0}.sh
   done
done

