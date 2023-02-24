#!/usr/bin/sh

now=`date +'%Y_%m_%d'`
mkdir -p LOGS/$now

if [ -e LOGS/$now/all_csvs ]; then
  echo LOGS/$now/all_csv exists
  else
  echo 'find all status.csv'
  find ./MONOMERS    -name 'status.csv' | grep -v -e TEST -e FAIL -e CHECK -e STOP >  LOGS/$now/all_csvs
  find ./DIMERS      -name 'status.csv' | grep -v -e TEST -e FAIL -e CHECK -e STOP >> LOGS/$now/all_csvs
  find ./CONVERGENCE -name 'status.csv' | grep -v -e TEST -e FAIL -e CHECK -e STOP >> LOGS/$now/all_csvs
  fi

grep -L "TERMINATED.NORMALLY,False," `cat LOGS/$now/all_csvs` > LOGS/$now/uncompleted_csv.dat

echo results in LOGS/$now/uncompleted_csv.dat
echo ">> wc results in LOGS/$now/uncompleted_csv.dat"
wc -l LOGS/$now/uncompleted_csv.dat
echo ">> echo results in LOGS/$now/uncompleted_csv.dat"
head LOGS/$now/uncompleted_csv.dat

echo 'These errors were found'
cat `cat LOGS/$now/uncompleted_csv.dat ` | sed 's/[0-9]/ /g; s/False/ /g; s/,.,/,/g; s/ //g; s/,,/,/g' | sort | uniq

grep 'exetyp=check recommended' `cat LOGS/$now/uncompleted_csv.dat ` | sed "s/status/ /g" | awk '{print $1}' > /data/mdi0316/WORK/IONIC_LIQUIDS/LOGS/$now/csv_error_check.dat
echo "`wc -l /data/mdi0316/WORK/IONIC_LIQUIDS/LOGS/$now/csv_error_check.dat` 'check recommended' written in /data/mdi0316/WORK/IONIC_LIQUIDS/LOGS/$now/csv_error_check.dat"

grep 'missing.input' `cat LOGS/$now/uncompleted_csv.dat ` | sed "s/status/ /g" | awk '{print $1}' > /data/mdi0316/WORK/IONIC_LIQUIDS/LOGS/$now/csv_missing_input.dat
echo "`wc -l /data/mdi0316/WORK/IONIC_LIQUIDS/LOGS/$now/csv_missing_input.dat` 'missing_input' written in /data/mdi0316/WORK/IONIC_LIQUIDS/LOGS/$now/csv_missing_input.dat"

grep 'missing.dir' `cat LOGS/$now/uncompleted_csv.dat ` | sed "s/status/ /g" | awk '{print $1}' > /data/mdi0316/WORK/IONIC_LIQUIDS/LOGS/$now/csv_missing_dir.dat
echo "`wc -l /data/mdi0316/WORK/IONIC_LIQUIDS/LOGS/$now/csv_missing_dir.dat` 'missing_dir' written in /data/mdi0316/WORK/IONIC_LIQUIDS/LOGS/$now/csv_missing_dir.dat"

grep 'dawrit' `cat LOGS/$now/uncompleted_csv.dat ` | sed "s/status/ /g" | awk '{print $1}' > /data/mdi0316/WORK/IONIC_LIQUIDS/LOGS/$now/csv_dawrit.dat
echo "`wc -l /data/mdi0316/WORK/IONIC_LIQUIDS/LOGS/$now/csv_dawrit.dat` 'dawrit' written in /data/mdi0316/WORK/IONIC_LIQUIDS/LOGS/$now/csv_dawrit.dat"

grep 'memory' `cat LOGS/$now/uncompleted_csv.dat ` | sed "s/status/ /g" | awk '{print $1}' > /data/mdi0316/WORK/IONIC_LIQUIDS/LOGS/$now/csv_memory.dat
echo "`wc -l /data/mdi0316/WORK/IONIC_LIQUIDS/LOGS/$now/csv_memory.dat` 'memory' written in /data/mdi0316/WORK/IONIC_LIQUIDS/LOGS/$now/csv_memory.dat"

