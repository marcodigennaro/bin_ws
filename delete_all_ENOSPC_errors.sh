#!/usr/bin/sh

folder=${1:-"/data/mdi0316/WORK/"}

echo 'Removing all ENOSPC folders in ' $folder

find $folder -name "log*" > all_log 
for ii in `cat all_log `;do grep -l ENOSPC $ii | sed "s/OPT/ /g" | awk '{print $1}' ;done > all_enospc_logs

echo `cat all_enospc_logs | wc -l` " found, remove? (Y)"
read -p "Proceed? (Y)" option
if [ $option == 'Y' ]; then
  rm -r `cat all_enospc_logs` 
fi

