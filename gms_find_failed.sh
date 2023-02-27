#!/usr/bin/sh

this_folder=`pwd`

if [ $# -eq 0 ]; then
   echo no argument provided
   folder=$this_folder
else
   folder=$this_folder/$1
   fi

echo "present folder = $this_folder"
echo "checking folders in : $folder"

rm -r `find  $folder -type d -name "GMS.FAIL"`

mkdir $folder/GMS.FAIL

find $folder -name '*inp'   | grep -v FAILED | grep -v CHECK | grep -v EQUIL | grep -v LOGS | grep -v NWCHEM | grep -v TEST > $folder/GMS.FAIL/all_inps
echo "`wc -l $folder/GMS.FAIL/all_inps | awk '{print $1}' ` inp files found in $folder"

find $folder -name 'log*'   | grep -v FAILED | grep -v CHECK | grep -v EQUIL | grep -v LOGS | grep -v NWCHEM | grep -v TEST > $folder/GMS.FAIL/all_logs
echo "`wc -l $folder/GMS.FAIL/all_logs | awk '{print $1}' ` log files found in $folder"

find $folder -name 'slurm*' | grep -v FAILED | grep -v CHECK | grep -v EQUIL | grep -v LOGS | grep -v NWCHEM | grep -v TEST > $folder/GMS.FAIL/all_slurm
echo "`wc -l $folder/GMS.FAIL/all_slurm | awk '{print $1}'` slurm files found in $folder"

for abs_log in `cat $folder/GMS.FAIL/all_logs`;  do

    #abs_log=$this_folder/$log

    if grep -q 'TERMINATED NORMALLY' $abs_log ; then
       echo $abs_log  >> $folder/GMS.FAIL/log_normally 
    else
       echo $abs_log  >> $folder/GMS.FAIL/log_not_normally 

       if   grep -q 'TERMINATED ABNORMALLY' $abs_log; then 
            echo $abs_log >> $folder/GMS.FAIL/log_abnormally
       elif grep -q 'ERROR: MEMORY REQUEST EXCEEDS AVAILABLE MEMORY' $abs_log; then
            echo $abs_log >> $folder/GMS.FAIL/log_not_normally_memory_error
       elif grep -q 'semget errno=ENOSPC -- check system limit for sysv semaphores.'  $abs_log; then
            echo $abs_log >> $folder/GMS.FAIL/log_not_normally_semaphore
       elif grep -q '0OPTIMIZATION ABORTED' $abs_log; then
            echo $abs_log >> $folder/GMS.FAIL/log_not_normally_opt_aborted
       elif grep -q 'FAILURE TO LOCATE STATIONARY POINT' $abs_log; then
            echo $abs_log >> $folder/GMS.FAIL/log_not_normally_opt_failed
       elif grep -q 'ddikick.x: Execution terminated due to error(s)'  $abs_log; then
            echo $abs_log >> $folder/GMS.FAIL/log_not_normally_error
       elif grep -q 'EINVAL'  $abs_log; then
            echo $abs_log >> $folder/GMS.FAIL/log_not_normally_einval
       elif grep -q 'Please save, rename, or erase these files from a previous run' $abs_log; then
            echo $abs_log >> $folder/GMS.FAIL/log_not_normally_resubmit
       elif grep -q 'EXETYP=CHECK is recommended' $abs_log; then
            echo $abs_log >> $folder/GMS.FAIL/log_not_normally_check
       else
            echo $abs_log >> $folder/GMS.FAIL/log_not_normally_unknown
       fi
    fi
    done

for slurm in `cat $folder/GMS.FAIL/all_slurm`;  do
    grep -l 'CANCELED'              $slurm >> $folder/GMS.FAIL/slurm_canceled
    grep -l 'slurmstepd: error:'    $slurm >> $folder/GMS.FAIL/slurm_stepd
    done

wc $folder/GMS.FAIL/*

for log_name in log_not_normally_semaphore log_not_normally_einval log_not_normally_unknown; do
    log_file="$folder/GMS.FAIL/$log_name"
    echo $log_file
    if [ -f "$log_file" ]; then
       wc $log_file
       read -p "Resubmit all? : " resubmit
       echo $resubmit
       if [ "$resubmit" == 'Y' ] || [ "$resubmit" == 'y' ]; then
          for ii in `cat $log_file`; do 
              loc_folder=`echo $ii | sed "s/log/ /g" | awk '{print $1}'`
              echo $loc_folder
              cd $loc_folder
              sh /home/mdi0316/bin/gms_submit_input.sh
              cd $this_folder
              done
          fi    
    rm $log_file
    else echo "does not exists"
       fi    
    echo ""
    done


log_file="$folder/GMS.FAIL/log_not_normally_memory_error"
echo $log_file
if [ -f "$log_file" ]; then
   wc $log_file
   read -p "Resubmit all? : " resubmit
   echo $resubmit
   if [ "$resubmit" == 'Y' ] || [ "$resubmit" == 'y' ]; then
      for ii in `cat $log_file`; do 
          loc_folder=`echo $ii | sed "s/log/ /g" | awk '{print $1}'`
          echo $loc_folder
          cd $loc_folder
          if [ -d "CHECK" ]; then
             echo CHECK folder exists in $loc_folder 
             echo $loc_folder >> $this_folder/tmp_check.dat
          else
             sh /home/mdi0316/bin/gms_make_check.sh
             fi
          cd $this_folder
          done
      fi    
fi
cp $this_folder/tmp_check.dat $log_file

#
#  log_num=`cat $folder/tmp_log | wc -l`
#  slurm_num=`cat $folder/tmp_slurm | wc -l`
#
#  if [ $log_num = 0 ]; then
#     echo $gms_dir has no logs 
#     echo $gms_dir has no logs >> $folder/find_no_complete/has_no_log 
#     inp=`ls $gms_dir | grep inp`
#     running=`qq | grep $inp | wc -l `
#     if [ $running = 0 ]; then
#        #cd $gms_dir
#        #cp /home/mdi0316/SCRIPTS/submit_gamess.sh ./
#        #sbatch -J $inp submit_gamess.sh $inp
#        #cd -
#        echo $gms_dir >> $folder/find_no_complete/no_log_no_running
#     elif [ $running = 1 ]; then
#        echo $gms_dir >> $folder/find_no_complete/running
#     else
#        echo "  WTF!!! $gms_dir "
#     fi
#
#  elif [ $log_num = 1 ]; then
# 
#     cat $folder/tmp_log
#
#     grep -l 'TERMINATED NORMALLY'   `cat $folder/tmp_log`   >> $folder/find_no_complete/terminated_normally
#     grep -l 'TERMINATED ABNORMALLY' `cat $folder/tmp_log`   >> $folder/find_no_complete/terminated_abnormally
#     grep -l 'semaphore'             `cat $folder/tmp_log`   >> $folder/find_no_complete/semaphore
#     grep -l 'EINVAL'                `cat $folder/tmp_log`   >> $folder/find_no_complete/einval
#     grep -l 'CANCELED'              `cat $folder/tmp_slurm` >> $folder/find_no_complete/canceled
#
#     opt=`echo $gms_dir | grep "/OPT/"`
#     if [ $? = 0 ]; then
#        grep -L 'LOCATED' `cat $folder/tmp_log` >> $folder/find_no_complete/opt_no_completed
#        fi
#
#  else
#     echo $gms_dir has two or more logs
#     echo $gms_dir >> $folder/find_no_complete/too_many_logs_dirs
#
#  fi
#
#  find $gms_dir -name 'status.csv' > $folder/tmp_csv
#  if [ $? != 0 ]; then 
#    grep ERROR `cat $folder/tmp_csv` >> $folder/find_no_complete/error_csv
#    fi
#
#done
