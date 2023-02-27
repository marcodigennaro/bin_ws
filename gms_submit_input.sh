#!/usr/bin/sh

inp_file=`ls *inp`

squeue | grep $inp_file

if [ $? = 0 ]; then
   echo $inp_file is running
   squeue | grep $inp_file
   else
   cp /home/mdi0316/SCRIPTS/submit_gamess.sh .
   echo "rm log* *dat *wfn slurm* "
   rm log* *dat *wfn slurm* 
   echo "sbatch -J $inp_file submit_gamess.sh $inp_file"
   sbatch -J $inp_file submit_gamess.sh $inp_file
   fi
