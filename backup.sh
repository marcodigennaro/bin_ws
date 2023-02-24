#!/bin/bash

month=`date +'%m_%b'`
MONTH=${month^^}

workstation_home="/home/mdi0316"

cluster_home="/home/mdi0316"
cluster_data="/data/mdi0316"
cluster_work="/data/mdi0316/WORK"
cluster_storage="/storage/mdi0316"

#for proj in MD-Batteries; do # QML IONIC_LIQUIDS MAGNETIC_MATERIALS  
#  sour=$cluster_work/$proj/src
#  cd $sour/
#  git add --all
#  git commit -m "auto_`date`"
#  git remote add origin https://cde.toyota-europe.com/stash/scm/~mdi0316/${proj}.git
#  git push -u origin master
#  cd

for fold in bin Inputfiles CLASSES FUNCTIONS SCRIPTS; do
  sour=$cluster_home/$fold
  for backup_dir in $cluster_home $cluster_data $cluster_storage; do
    dest=$backup_dir/BACKUP/CLUSTER/HOME/$fold
    echo backing $sour/ to $dest
    mkdir -p $dest
    rsync -ra --update $sour/ $dest
    done
  rsync -ra --update $sour/ mdi0316@10.100.192.47:/home/mdi0316/BACKUP/CLUSTER/HOME/$fold
done

for proj in MWCNT MD-Batteries QML IONIC_LIQUIDS MAGNETIC_MATERIALS; do
  for tmp_dir in src data CSV JSON LOGS ; do
    sour=$cluster_work/$proj/$tmp_dir
    if [ -e $sour ]; then
      for backup_dir in $cluster_home $cluster_data $cluster_storage; do
        dest=$backup_dir/BACKUP/CLUSTER/WORK/$proj/$tmp_dir
        echo backing $sour/ to $dest
        mkdir -p $dest
        rsync -ra --update $sour/ $dest
        done
      rsync -ra --update $sour/ mdi0316@10.100.192.47:/home/mdi0316/BACKUP/CLUSTER/WORK/$proj/$tmp_dir
    fi
  done

  ## rsync RUNS only on /storage
  #sour=$cluster_work/$proj/RUNS
  #dest=$cluster_storage/BACKUP/CLUSTER/WORK/$proj/RUNS
  #echo backing $sour/ to $dest
  #mkdir -p $dest
  #rsync -ra --update $sour/ $dest

  done
