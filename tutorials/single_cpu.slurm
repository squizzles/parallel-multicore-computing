#/bin/bash /Users/rob/Documents/University/COMP90020/distributed-mutual-exclusion/run.sh

#SBATCH --job-name=singlecpu
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1

# Run two tasks (separate CPUs)
srun --ntasks=1 echo "Task 1"
srun --ntasks=1 echo "Task 2"

