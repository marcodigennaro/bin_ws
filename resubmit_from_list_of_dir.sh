#!/usr/bin/sh

list_of_dirs=$1

for dir in `cat $1`; do
  cd $dir
  inp=`ls | grep inp`
  sbatch -J $inp -p nodeshiq submit_gamess.sh $inp
  cd -
  done
