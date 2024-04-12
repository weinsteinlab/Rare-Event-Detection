#!/bin/bash

## declare an array variable
StoragePath='/athena/hwlab/midtier/hex4001/wasabi/StarD4-chol-phase2-batch2/raw_swarms/' #0..35,49
system_name='CHLp2b2'
further_stride=5
scriptfolder='/home/hex4001/RED'
## now loop through the above array
for replica in {0..35}
#for replica in 6
do
   replica0=`printf %04d $replica`
   for j in {0..49}
   #for j in 0
   do
      frame=`expr $j \* 1000`
      frame0=`printf %06d $frame`
      aimname="data/${system_name}_traj${replica0}/${system_name}_traj${replica0}_${frame0}.gz"
      if test -f "$aimname"; then
         echo "$system_name traj$replica0 frame$frame0 done"
         #continue
      fi
      sbatch ${scriptfolder}/submit_contact_from_storage.sh $StoragePath $system_name $replica0 $j $further_stride
   done
done

