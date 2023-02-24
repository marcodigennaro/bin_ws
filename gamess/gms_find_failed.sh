#!/usr/bin/sh

this_folder=`pwd`
now=`date +'%Y_%m_%d'`

find_files()
{
  echo "find $1 -name $2 | grep -v -e PCseg-4 -e CC6 -e CC5 -e MOFS -e FAILED -e CHECK -e STOPPED -e LOGS -e TEST > $3 "
        find $1 -name $2 | grep -v -e PCseg-4 -e CC6 -e CC5 -e MOFS -e FAILED -e CHECK -e STOPPED -e LOGS -e TEST > $3
}

get_dir_from_abs_log()
{
  echo $1 | sed "s/log/ /g" | awk '{print $1}'
}

if [ $# -eq 0 ]; then
   echo no argument provided
   folder=$this_folder
else
   folder=$this_folder/$1
   fi

echo "present folder = $this_folder"
echo "checking folders in : $folder"

echo $now
#exit
#rm -r `find  $folder -type d -name "GMS.FAIL/$now"`

now_dir="/data/mdi0316/WORK/IONIC_LIQUIDS/LOGS/$now"
log_dir="/data/mdi0316/WORK/IONIC_LIQUIDS/LOGS/$now/GMS.FAIL"
mkdir -p $log_dir

if [ -e $log_dir/all_inps ]; then
   echo $log_dir/all_inps exists 
   else
   find_files $folder '*inp' $log_dir/all_inps
   fi

if [ -e $log_dir/all_logs ]; then
   echo $log_dir/all_logs exists 
   else
   find_files $folder 'log*' $log_dir/all_logs
   fi

if [ -e $log_dir/all_slurms ]; then
   echo $log_dir/all_slurms exists 
   else
   find_files $folder 'slurm*' $log_dir/all_slurms
   fi

if [ -e $log_dir/log_normally ]; then
   echo $log_dir/log_normally exists 
   else
   echo "grep -l 'TERMINATED NORMALLY' \`cat $log_dir/all_logs\` > $log_dir/log_normally "
         grep -l 'TERMINATED NORMALLY' `cat $log_dir/all_logs` > $log_dir/log_normally 
   fi

if [ -e $log_dir/log_abnormally ]; then
   echo $log_dir/log_abnormally exists 
   else
   echo "grep -l 'TERMINATED -ABNORMALLY-' \`cat $log_dir/all_logs\` > $log_dir/log_abnormally "
         grep -l 'TERMINATED -ABNORMALLY-' `cat $log_dir/all_logs` > $log_dir/log_abnormally 
   fi

echo "`wc -l $log_dir/all_inps | awk '{print $1}' ` inp files found in $folder"
echo "`wc -l $log_dir/all_logs | awk '{print $1}' ` log files found in $folder"
echo "`wc -l $log_dir/all_slurms | awk '{print $1}'` slurm files found in $folder"
echo "`wc -l $log_dir/log_normally | awk '{print $1}'` terminated normally files found in $folder"
echo "`wc -l $log_dir/log_abnormally | awk '{print $1}'` terminated not normally files found in $folder"

rm $log_dir/log_not_normally_einval
rm $log_dir/log_not_normally_error
rm $log_dir/log_not_normally_memory
rm $log_dir/log_not_normally_gradient
rm $log_dir/log_not_normally_dawrit
rm $log_dir/log_not_normally_principal_axis
rm $log_dir/log_not_normally_remapping
rm $log_dir/log_not_normally_unknown

for abs_log in `cat $log_dir/log_abnormally`;  do

    # SEMAPHORES/EINVAL ERROR
    if grep -q -e 'check system limit for sysv semaphores' \
               -e 'EINVAL'  $abs_log; then
         echo $abs_log >> $log_dir/log_not_normally_einval

    # GAMESS ERROR
    elif grep -q -e '0OPTIMIZATION ABORTED' \
                 -e 'FAILURE TO LOCATE STATIONARY POINT' \
                 -e ' SCF DID NOT CONVERGE...NO MPLEVL=2 CALCULATION' \
                 -e 'SCF IS UNCONVERGED, TOO MANY ITERATIONS' \
                 -e 'Please save, rename, or erase these files from a previous run' $abs_log; then
         echo $abs_log >> $log_dir/log_not_normally_error

    # MEMORY ERROR
    elif grep -q -e 'INSUFFICIENT REPLICATED MEMORY REQUESTED' \
                 -e 'EXETYP=CHECK is recommended' \
                 -e 'ERROR: MEMORY REQUEST EXCEEDS AVAILABLE MEMORY' $abs_log; then
         echo $abs_log >> $log_dir/log_not_normally_memory
         get_dir_from_abs_log $abs_log >> $now_dir/memory_check.dat

    elif grep -q 'NUMERICAL GRADIENT CANNOT PROCEED' $abs_log; then
         echo $abs_log >> $log_dir/log_not_normally_gradient

    elif grep -q 'DAWRIT' $abs_log; then
         echo $abs_log >> $log_dir/log_not_normally_dawrit

    elif grep -q 'UNABLE TO GENERATE PRINCIPAL AXES' $abs_log; then
         echo $abs_log >> $log_dir/log_not_normally_principal_axis

    elif grep -q 'I SHELLS FROM GAMESS ORDER' $abs_log; then
         echo $abs_log >> $log_dir/log_not_normally_remapping

    else
         echo $abs_log >> $log_dir/log_not_normally_unknown
    fi
    done

#for slurm in `cat $log_dir/all_slurm`;  do
#    grep -l 'CANCELED'              $slurm >> $log_dir/slurm_canceled
#    grep -l 'slurmstepd: error:'    $slurm >> $log_dir/slurm_stepd
#    done
#
#wc $log_dir/*
#
#for log_name in log_not_normally_semaphore log_not_normally_einval log_not_normally_unknown; do
#    log_file="$log_dir/$log_name"
#    echo $log_file
#    if [ -f "$log_file" ]; then
#       wc $log_file
#       read -p "Resubmit all? : " resubmit
#       echo $resubmit
#       if [ "$resubmit" == 'Y' ] || [ "$resubmit" == 'y' ]; then
#          for ii in `cat $log_file`; do #              loc_folder=`echo $ii | sed "s/log/ /g" | awk '{print $1}'`
#              loc_folder=`echo $ii | sed "s/log/ /g" | awk '{print $1}'`
#              echo $loc_folder
#              cd $loc_folder
#              sh /home/mdi0316/bin/gms_submit_input.sh
#              cd $this_folder
#              done
#          fi    
#    rm $log_file
#    else echo "does not exists"
#       fi    
#    echo ""
#    done

find $log_dir -type f -empty -delete
