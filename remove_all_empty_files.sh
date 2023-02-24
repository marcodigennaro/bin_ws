#!/bin/bash

if [ -z "$1" ]; then
  fold="./"
else
  fold="$1"
fi

while true; do

  echo specified folder = $fold
  echo specified depth "(\$2)" = $2
  if [ -z "$2" ]; then
    echo "Finding empty files in $fold"
    all_empty=`find $fold -type f -empty`
    find $fold -type f -empty > all_empty.txt
  else
    echo "Finding empty files in $fold at depth $2"
    all_empty=`find $fold -mindepth $2 -maxdepth $2 -type f -empty`
    find $fold -mindepth $2 -maxdepth $2 -type f -empty > all_empty.txt
  fi

  [ -s "all_empty.txt" ]
  if [ $? == 0 ]; then
    while IFS= read -r line
    do
      echo '  Removing ' $line
      rm -r "$line"
    done < "all_empty.txt"
  else
    echo '  No empty files was found'
    break
  fi

  done

