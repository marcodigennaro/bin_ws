#!/usr/bin/sh

now=`date +'%Y_%m_%d'`
mkdir -p LOGS/$now/CSV

for mono in PF6 BF4 EMIM; do
  echo "Running mono << $mono >>"
  ./monomers.py $mono | tee LOGS/$now/log_$mono
  cp MONOMERS/CSV/${mono}*csv LOGS/$now
  done

for dime in EMIM_BF4 EMIM_PF6; do
  echo "Running dime << $dime >>"
  ./dimers.py $dime | tee LOGS/$now/log_$dime
  cp -r DIMERS/${dime}/CSV/ LOGS/$now/CSV/$dime
  done

sh ~/bin/gms_find_all_memory_error_logs.sh  
sh ~/bin/gms_find_all_uncompleted_csv.sh

