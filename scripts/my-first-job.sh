#!/bin/bash
#SBATCH -t 10:00
#SBATCH -c 1
#SBATCH --constraint=skylake

echo "================================================================="
echo "                        Welcome to Spider"
echo "================================================================="
echo ""
echo "The status of the worker nodes is as displayed below:"
sinfo
echo ""
echo "Your current running jobs are listed below:"
squeue -u $USER
echo ""
echo "You just ran your first job on" $HOSTNAME " with a job ID " $SLURM_JOBID
echo ""
sleep 15s
