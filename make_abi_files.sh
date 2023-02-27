#!/bin/bash

ha2ev=27.2114
main_fold="/home/mdi0316/WORK/SEP/ABINIT/EMIM_BF4/"

bf4_log="/home/mdi0316/WORK/SEP/ABINIT/BF4/log"
bf4_en_ha=`grep "etotal" $bf4_log | awk '{printf "%4.4f", $2}'`
bf4_en_ev=$(expr $bf4_en_ha*$ha2ev | bc )

emim_log="/home/mdi0316/WORK/SEP/ABINIT/EMIM/log"
emim_en_ha=`grep "etotal" $emim_log | awk '{printf "%4.4f", $2}'`
emim_en_ev=$(expr $emim_en_ha*$ha2ev | bc )

for Z in `seq 10 19`; do
  tot_en_ev='nan'
  int_en_ev='nan'
  new_fold="$main_fold/SCAN/Z_$Z"
  cd $new_fold
  calc_compl=`grep " Calculation completed." log_$Z | wc -l`
  if [ $calc_compl == 1 ]; then
    tot_en_ha=`grep "etotal" log_$Z | awk '{printf "%4.4f", $2}'`
    tot_en_ev=$(expr $tot_en_ha*$ha2ev | bc )
    int_en_ev=`echo $tot_en_ev  $emim_en_ev $bf4_en_ev | awk '{print $1-$2-$3}' `
    fi
  echo $Z $int_en_ev
  #cp emim_bf4.in emim_bf4.files $new_fold
  #abinit < emim_bf4.files | tee log_$Z
  cd $main_fold
  done

