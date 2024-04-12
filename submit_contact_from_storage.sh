#!/bin/bash
#SBATCH -p edison,panda_physbio
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --exclude=node150

vmd3=/home/des2037/vmd3/bin/vmd3
scriptfolder='/home/hex4001/RED'
StoragePath=$1
system_name=$2
trajid=$3
startslice=$4
further_stride=$5

echo "inputs: $StoragePath $system_name $trajid $startslice $further_stride"

startf=`expr $startslice \* 1000`
section_length=1000
endf=`expr ${startf} + ${section_length} - 1`
startf0=`printf %06d $startf`

psfname=`ls $StoragePath/swarm0000/swarm0000_traj${trajid}/*psf`
xtcname=`ls $StoragePath/swarm0000/swarm0000_traj${trajid}/*xtc`
logname=`find $StoragePath/swarm0000/swarm0000_traj${trajid}/*log -name "*subjob*"`

#find out the subjobs files that contribute to the "startf"-to-"endf" timezone
framecount=0
reading=0
readstartf=Null
readfilelist=''
for i in $logname
do
  fileframe=`wc ${i%.???}.log -l | awk '{printf "%d ",$1}'`
  if [ "$fileframe" -le 1 ] ; then
    continue
  fi
  fileframe=`expr $fileframe - 1`
  framecount=`expr $framecount + $fileframe`
  if [ "$framecount" -gt "$startf" ] ; then
    if [ "$reading" -eq 0 ] ; then
      reading=1
      readstartf=`expr $fileframe - $framecount + $startf`
    fi
    readfilelist="$readfilelist ${i%.???}.xtc"
    if [ $framecount -gt $endf ] ; then
      break
    fi
  fi
done
readendf=`expr $readstartf + $section_length - 1`


mkdir -p data/${system_name}_traj${trajid}
mkdir check

echo "$psfname"
echo "$readfilelist"
echo "$readstartf data/${system_name}_traj${trajid}/${system_name}_traj${trajid}_${startf0}.dat $further_stride"

$vmd3 -dispdev text $psfname $readfilelist -e ${scriptfolder}/first_contact_from_storage.tcl -args $readstartf data/${system_name}_traj${trajid}/${system_name}_traj${trajid}_${startf0}.dat $further_stride
python3 ${scriptfolder}/process_contacts.py data/${system_name}_traj${trajid}/${system_name}_traj${trajid}_${startf0}.dat data/${system_name}_traj${trajid}/${system_name}_traj${trajid}_${startf0}
rm data/${system_name}_traj${trajid}/${system_name}_traj${trajid}_${startf0}.dat

exit
