#!/bin/bash -l

#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=nodesloq
#SBATCH --job-name=find_directories
#SBATCH --output=out.%j
#SBATCH --error=err.%j

start_date=`date`
start_time="$(date -u +%s.%N)"
echo "Job slurm id: $SLURM_JOB_ID"
echo "Running on:   $SLURM_NODELIST"
echo "    nnodes       - $SLURM_NNODES"
echo "    processors   - $SLURM_NPROCS"
echo "    ntasks       - $SLURM_NTASKS"
echo "    cpu-per-task - $SLURM_CPUS_PER_TASK"

# load modules
module purge
module load slurm

find /data/mdi0316 -type d > ~/all_directories.txt
touch ~/all_sizes.txt
truncate -s0 ~/all_sizes.txt

while IFS="" read -r p || [ -n "$p" ]
do
  find $p -type f -maxdepth 1 > ~/tmp_file.txt
  size=`wc -l tmp_file.txt`
  echo $size $p  >> ~/all_sizes.txt
done < ~/all_directories.txt

end_date=`date`
end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"
echo "Starting at: $start_date"
echo "Ending at:   $end_date"

echo "Total of $elapsed seconds elapsed for process"


# CommonAdapter (SLURM) completed writing Template
