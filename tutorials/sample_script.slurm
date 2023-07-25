#!/bin/bash

#SBATCH --job-name=sample-script
# maximum duration of the job, format hh:mm:ss
#SBATCH --time=00:00:60
# no other jobs using hardware
#SBATCH --exclusive
# total number of tasks (processes)
#SBATCH --ntasks=2
# number of nodes
#SBATCH --nodes=2
# multi-threaded -> 1 task per node; process-parallelism -> number cores per node
#SBATCH --tasks-per-node=1
# threads per process
#SBATCH --cpus-per-task=1

# parition := nodes job will execute on
# Spartan "physical" partition is for running on CPU nodes
#SBATCH --partition=physical
#SBATCH --qos=standard

# Run two tasks (separate CPUs)
srun --ntasks=1 echo "Task 1"
srun --ntasks=1 echo "Task 2"

