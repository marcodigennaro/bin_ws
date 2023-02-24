#!/usr/bin/sh

inp=`ls *inp`

sbatch -J $inp submit_gamess.sh $inp 
