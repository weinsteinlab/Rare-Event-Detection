#!/bin/bash

## declare an array variable
folders='
CHLp2b2
'
starttraj=0
endtraj=11
startslice=0
endslice=49
scriptfolder='/home/hex4001/RED'
## now loop through the above array
for i in $folders; do
   for replica in `seq $starttraj $endtraj`; do
      replica0=`printf %04d $replica`
      echo "$i $replica"
      cp $scriptfolder/submit_dataconcat.sh temp_${i}_${replica0}.sh
      sed -i "s/input3/$i/g" temp_${i}_${replica0}.sh
      sed -i "s/input4/$replica0/g" temp_${i}_${replica0}.sh
      sed -i "s/input1/$startslice/g" temp_${i}_${replica0}.sh
      sed -i "s/input2/$endslice/g" temp_${i}_${replica0}.sh
      sbatch temp_${i}_${replica0}.sh
      rm temp_${i}_${replica0}.sh
   done
done

