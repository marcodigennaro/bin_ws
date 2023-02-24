#!/usr/bin/sh

backup_dirs='
/home/mdi0316/BACKUP/CLUSTER/
/data/mdi0316/BACKUP/CLUSTER/
/storage/mdi0316/BACKUP/CLUSTER/ 
'
echo `date +'%H:%M - %m/%d/%Y'`
echo `date +'%H:%M - %m/%d/%Y'` > ~/.backup_files


## CASE 1: copy all files in source dir
source_dirs='
/home/mdi0316/bin
/home/mdi0316/SCRIPTS
/home/mdi0316/CLASSES
/home/mdi0316/FUNCTIONS
/home/mdi0316/Inputfiles
'

for bb in `echo $backup_dirs`; do 
  mkdir -p $bb/bin
  mkdir -p $bb/CLASSES
  mkdir -p $bb/FUNCTIONS
  mkdir -p $bb/Inputfiles
  mkdir -p $bb/WORK
  mkdir -p $bb/SCRIPTS
  done

for directory in `echo $source_dirs`; do
  for bb in `echo $backup_dirs`; do 
    echo rsyncing full directory $directory to $bb
    rsync --update -raz $directory $bb
    done
 done

## CASE 2: copy only .sh and .py files from source dir
source_dirs='
/data/mdi0316/WORK/
'
#/data/mdi0316/WORK/MOFS/
#/data/mdi0316/WORK/IONIC_LIQUIDS/
#/data/mdi0316/WORK/GAS_EXHAUST/
#/data/mdi0316/WORK/NWCHEM/
#/data/mdi0316/WORK/QC/

for directory in `echo $source_dirs`; do
    echo $directory
    for extension in py sh; do
        echo "find $directory  -maxdepth 5 -type f -name '*.$extension' | grep -v -e Activity_Test -e runorca -e submit -e decorated -e slurm -e quantum-machine-9 >> ~/.backup_files"
        find      "$directory" -maxdepth 5 -type f -name "*.$extension" | grep -v -e Activity_Test -e runorca -e submit -e decorated -e slurm -e quantum-machine-9 >> ~/.backup_files
        done
    done

for ff in `cat ~/.backup_files`; do 
  for bb in `echo $backup_dirs`; do 
    echo rsyncing file $ff to $bb/WORK/
    rsync --update -raz $ff $bb/WORK/
    done
 done

echo log in ~/.backup_files
