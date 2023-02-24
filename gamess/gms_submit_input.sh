#!/usr/bin/sh

inp_file=`ls *inp`

squeue | grep $inp_file

if [ $? = 0 ]; then
   echo $inp_file is running
   squeue | grep $inp_file
else
   if [ ! -e submit_gamess.sh ]; then
      cp /home/mdi0316/SCRIPTS/submit_gamess.sh .
   fi 
   echo "rm log* *dat *wfn slurm* *csv err.gms"
   rm log* *dat *wfn slurm* *csv err.gms
   echo "sbatch -J $inp_file submit_gamess.sh $inp_file"
   sbatch -J $inp_file submit_gamess.sh $inp_file
fi
