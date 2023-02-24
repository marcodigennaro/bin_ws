#!/usr/bin/sh

now=`date +'%Y_%m_%d'`
mkdir -p LOGS/$now

for mono in PF6 BF4 EMIM; do
  ./monomers.py $mono | tee LOGS/$now/log_$mono
  cp MONOMERS/CSV/${mono}*csv LOGS/$now
  done

#for mono in SCN C1MIM C3MIM C4MIM C5MIM C6MIM C7MIM C8MIM; do
#  ./monomers.py $mono | tee LOGS/log_$mono
#  done
