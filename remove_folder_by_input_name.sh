#!/usr/bin/sh

#$1 can be name or slurm_ID

echo '============= START ============='
job=`squeue -u \`whoami\` -o "%.8i %.9P %.90j %.8T %.10M %.9l %.6D %R" | grep $1 | awk '{print $1}'`
name=$1 
dir=`find $2 -name $1`

echo "job=$job"
echo "name=$name"
echo "dir=$dir"

if   [ "$dir" = '' ]; then
   echo folder not found
elif [ "$job" = '' ]; then
   echo jobid not found
else
   echo folder found
   dir=`find $2 -name $1 | sed "s/$1//g"`
   echo $dir
   echo DELETING JOB:$job NAME: $1
   scancel $job
   echo REMOVING $dir
   rm -rf $dir
fi
echo '-------------- END --------------'
