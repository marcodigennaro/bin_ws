#!/usr/bin/sh

squeue -u mdi0316 | awk '{print $1}' > ~/.tmp_running_idx

all_ids=`ls /data/scratch-no-backup/mdi0316/rungms/`

for idx in $all_ids; do 
  running=`cat ~/.tmp_running_idx | grep $idx`
  if [ $? = 0 ]; then
     echo $idx is running 
  else
     echo $idx is not running 
     rm -r /data/scratch-no-backup/mdi0316/rungms/$idx
     fi
  done
