#!/bin/bash

bat=task.bat
#run=/home/mdi0316/Inputfiles/EASY/NKK_run_lammps.py
dat=`ls *task.dat`
lmp=`ls *lammps.dat`

echo $bat
#echo $run
echo $dat
echo $lmp

rm -r easy_quick_test
mkdir easy_quick_test

cp ${bat} easy_quick_test/
#cp ${run} easy_quick_test/
cp ${lmp} easy_quick_test/
cp ${dat} easy_quick_test/task.dat

cd easy_quick_test

./${bat}

cat task.res
cat delta.csv
