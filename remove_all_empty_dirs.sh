#!/usr/bin/sh

if [ -z "$1" ]; then
  fold="./"
else
  fold="$1"
fi

while true; do
  echo "Finding empty folders in $fold"
  all_empty=`find $fold -type d -empty`
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

