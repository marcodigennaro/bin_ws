#!/usr/bin/bash


for ii in TRIM NR RFO QA SCHLEGEL CONOPT;do mkdir $ii ; cd $ii ; sed "s/NSTEP=1000/NSTEP=1000\n  METHOD=$ii/g" ../default.inp > $ii.inp ; cp ~/scripts/submit_gamess.sh . ; sbatch -p nodeshiq -J $ii submit_gamess.sh $ii;  cd - ;done

