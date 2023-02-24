#!//bin/sh


count=1
cat $1 | sort | uniq > tmp_list

tot=`cat tmp_list | wc -l`

for write_dir in `cat tmp_list`;do
  echo "===================== $count/$tot"
  echo Fixing $write_dir
  inp_file=`ls $write_dir | grep inp`
  echo $inp_file
  qq | grep $inp_file
  if [ $? == 0 ]; then
      echo '>> Running, skipping'
  else
      p_dir=` echo $write_dir | sed "s/\// /g" | awk '{$NF=""; print $0}' | sed 's/ /\//g' `
      echo P_dir: $p_dir
      grep 'LOCATED' /$p_dir/R_*/log* | sort
      echo Fixing $write_dir
      ls $write_dir/OPT.STOPPED
      read -p '>> choose R to read from : ' copy_R
      read_log=`ls /$p_dir/R_${copy_R}/log* `
      echo '>> will read from: ' $read_log
      grep LOCATED $read_log
      read -p '>> proceed? : ' proceed
      if [ $proceed == 'Y' ]; then
        cd $write_dir
        python ~/bin/gms_resubmit_optimization.py $read_log
        cd /data/mdi0316/WORK/IONIC_LIQUIDS
        fi
  fi
  count=$(( $count + 1 ))
  echo ''
done
