#!/usr/bin/sh


count=1
cat $1 | sort | uniq > tmp_list

tot=`cat tmp_list | wc -l`

for ii in `cat tmp_list`;do cd $ii; echo ==== $count/$tot === ; gms_make_check.py ; count=$(( count + 1 )); cd /data/mdi0316/WORK/IONIC_LIQUIDS/ ; echo '' ; done

