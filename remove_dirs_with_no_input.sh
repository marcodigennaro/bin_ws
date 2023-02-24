#!/usr/bin/sh

if [ -z "$1" ]; then
  fold="./"
else
  fold="$1"
fi

rm dirs_with_no_input
rm dirs_with_subfolders
echo looking in $fold
find $fold -type d -name "R_*" | grep SCAN | grep RUNS | grep -v SCAN_from_ISOLATED > all_R_dirs
echo `wc -l all_R_dirs` directory found

for ii in `cat all_R_dirs`;do 
  ls $ii | grep inp > /dev/null
  if [ $? != 0 ];  then
   echo $ii >> dirs_with_no_input
   fi
  done 


# all folders with subfolders
for ii in `cat dirs_with_no_input`;do find $ii -mindepth 1 -type d ;done >> dirs_with_subfolders


echo `wc -l dirs_with_no_input`   dirs with no input found
echo `wc -l dirs_with_subfolders` of these have subfolders

