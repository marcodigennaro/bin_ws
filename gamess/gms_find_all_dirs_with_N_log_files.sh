#!/usr/bin/sh

this_folder=`pwd`
now=`date +'%Y_%m_%d'`
if [ $# -eq 0 ]; then
   echo no argument provided
   folder=$this_folder
else
   folder=$this_folder/$1
   fi

echo "present folder = $this_folder"
echo "checking folders in : $folder"

echo $now

all_basis='STO N21 N31 N311 DZV CCD CCT CCQ G3L TZV PCseg-0 PCseg-1 PCseg-2 PCseg-3 APCseg-0 APCseg-1 APCseg-2 APCseg-3'

find $folder -type d | grep -v -e MOFS -e CHECK -e STOP -e FAIL -e TEST -e LOG > all_dirs
for dir in `cat all_dirs`; do
  last_dir=`echo $dir | sed "s/\// /g" | awk '{print $NF}'`
  echo $all_basis | grep $last_dir 
  if [ $? != 0 ]; then
    num_log=`ls $dir/log* | wc -l`
    num_inp=`ls $dir/*inp | wc -l`
    if [[ $num_inp == 0 ]]; then
      echo no inp file $dir >> /data/mdi0316/WORK/IONIC_LIQUIDS/LOGS/$now/wrong_inp_log_number.dat
    elif [[ $num_log != 0 && $num_log != 1 ]]; then
      echo 2 or more logs $dir >> /data/mdi0316/WORK/IONIC_LIQUIDS/LOGS/$now/wrong_inp_log_number.dat
    elif [[ $num_log != 1 || $num_inp != 1 ]]; then
      echo $num_inp $num_log $dir >> /data/mdi0316/WORK/IONIC_LIQUIDS/LOGS/$now/wrong_inp_log_number.dat
      fi
    fi
  done

