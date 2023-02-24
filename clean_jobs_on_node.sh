#!/usr/bin/sh

while true ; do 
  for node in $1 ; do 
    id=`qq | grep $node | awk '{print $3}'` 
    if [ "$id" != '' ]; then 
       remove_folder_by_input_name.sh $id $2
    else:
       echo no running id
    fi  
    done
  sleep 10
  done
