#/usr/bin/sh

# usage
# delete_jobs.sh str


if [ $1 == 'all' ]; then
  list_of_jobs=`squeue -u \`whoami\` -o "%i %P %j %u %T" | awk '{print $1}'`
else
  list_of_jobs=`squeue -u \`whoami\` -o "%i %P %j %u %T" | grep $1 | awk '{print $1}'`
fi

length=`echo $list_of_jobs | wc `
echo "list of jobs:"
echo "$list_of_jobs"
echo "$length total jobs" 


read -p "Proceed? (Y)" option

if [ $option == 'Y' ]; then
  for job in $list_of_jobs; do
    echo "deleting $job"
    scancel $job
    done
  fi
