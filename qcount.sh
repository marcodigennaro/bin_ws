#/usr/bin/sh


count_running () {
  run=`squeue -u \`whoami\` -o "%.8i %.9P %.50j %.8u %.8T %.10M %.9l %.6D %R" | grep $1  | wc -l`
  echo "$run $1 jobs"
  return $run
}


squeue -u `whoami` -o "%.8i %.9P %.50j %.8u %.8T %.10M %.9l %.6D %R"

count_running "RUNNING"
count_running "PENDING"
count_running "nodeshiq"
count_running "nodesloq"

