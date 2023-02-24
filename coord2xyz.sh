#!/bin/bash

input=$1
while IFS= read -r line
do
  echo $line | grep '\$'
  if ! [ $? -eq 0 ]; then 
    echo "$line" | awk '{print $4"  "$1"  "$2"  "$3}'
  fi
done < "$input"
