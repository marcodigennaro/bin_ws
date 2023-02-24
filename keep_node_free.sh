#!/usr/bin/sh

work_dir=/data/mdi0316/WORK

node=node${1}
fold=${2}

inp=`qq | grep $node | awk '{print $3}'`

echo $inp | grep gms_ 

if [ $? == 0 ]; then

idx=`qq | grep $node | awk '{print $1}'`
echo "find $work_dir/$fold -name \"log*${idx}\""
log=`find $work_dir/$fold -name "log*${idx}"`
dir=`echo $log | sed "s/log/ /g" | awk '{print $1}' `

cd $dir
scancel $idx
gms_submit_input.sh
cd $work_dir
fi
