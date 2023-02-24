#!/usr/bin/sh

now=`date +'%Y_%m_%d'`
mkdir -p LOGS/$now

memory_check_file=/data/mdi0316/WORK/IONIC_LIQUIDS/LOGS/$now/memory_check.dat

push_folder()
{
    fold=`echo $1 | sed "s/log/ /g" | awk '{print $1}'`
    echo "$fold >> $memory_check_file"
    echo  $fold >> $memory_check_file
}

#if [ -e LOGS/$now/all_logs ]; then
#  echo LOGS/$now/all_logs exists
#  else
  echo 'find all logs'
  find ./MONOMERS    -name 'log.gms*' | grep -v TEST | grep -v FAIL | grep -v CHECK | grep -v old | grep -v test | grep -v STOP > LOGS/$now/all_logs
  find ./DIMERS      -name 'log.gms*' | grep -v TEST | grep -v FAIL | grep -v CHECK | grep -v old | grep -v test | grep -v STOP >> LOGS/$now/all_logs
  find ./CONVERGENCE -name 'log.cc*'  | grep -v TEST | grep -v FAIL | grep -v CHECK | grep -v old | grep -v test | grep -v STOP >> LOGS/$now/all_logs
#  fi

#echo `wc -l LOGS/$now/all_logs` log files found

#if [ -e LOGS/$now/all_abnormally_logs ]; then
#  echo LOGS/$now/all_abnormally_logs exists
#  else
  echo 'find all abnormally' 
  grep -l 'ABNORMALLY' `cat LOGS/$now/all_logs` > LOGS/$now/all_abnormally_logs
#  fi

echo `wc -l LOGS/$now/all_abnormally_logs` abnormally finished log files found

#if [ -e $memory_check_file ]; then
#  echo $memory_check_file exists
#  else
  echo "find all memory known message > $memory_check_file"
  for log in `cat LOGS/$now/all_abnormally_logs`; do
    grep 'ERROR: MEMORY REQUEST EXCEEDS AVAILABLE MEMORY'          $log
      if [ $? == 0 ]; then push_folder $log; fi
    grep 'THE DISTRIBUTED MEMORY REQUIRED FOR THIS STEP IS MEMDDI' $log  
      if [ $? == 0 ]; then push_folder $log; fi
    grep 'INSUFFICIENT DISTRIBUTED MEMORY REQUESTED'               $log
      if [ $? == 0 ]; then push_folder $log; fi 
    grep 'INSUFFICIENT DISTRIBUTED MEMORY REQUESTED'               $log
      if [ $? == 0 ]; then push_folder $log; fi
    grep 'PLEASE INCREASE -MWORDS- IN $SYSTEM APPROPRIATELY'       $log
      if [ $? == 0 ]; then push_folder $log; fi
    #grep -l 'ADD' $log | grep -l 'MORE WORDS'
    #grep -l 'MEMDDI =' $log | grep -l 'BUT NEEDS TO BE ='
    #grep -l 'PROCESS NO.' $log | grep -l 'REQUIRED='  
    done
#  fi

