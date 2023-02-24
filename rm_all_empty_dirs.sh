#!/bin/sh

while true; do
  echo 'Finding empty folders'
  all_empty=`find . -type d -empty`
  if [ -z "$all_empty" ]; then
    echo '  No empty folder was found'
    break
  else
    for ii in $all_empty; do
      echo '  Removing ' $ii
      rm -r $ii
      done
  fi
  done

