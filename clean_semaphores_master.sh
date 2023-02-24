#!/usr/bin/sh

pdsh -w node0[01-16] /home/mdi0316/bin/clean_semaphores_nodes.sh

echo "which folder?"

folder="$1"

echo 'Removing all ENOSPC/EINVAL folders in ' $folder

find $folder -name "log*" > all_log 

for ii in `cat all_log `;do grep -l ENOSPC $ii | sed "s/log/ /g" | awk '{print $1}' ;done > all_enospc_logs
for ii in `cat all_log `;do grep -l EINVAL $ii | sed "s/log/ /g" | awk '{print $1}' ;done > all_einval_logs


pdsh -w nodes0[01-15] /home/mdi0316/bin/clean_semaphores_nodes.sh

echo `cat all_enospc_logs all_einval_logs | wc -l` " found"

read -p "Proceed? (Y)" option
if [ $option == 'Y' ]; then
for file in `cat all_enospc_logs all_einval_logs`; do
  echo $file
  rm -r $file 
  done
fi
