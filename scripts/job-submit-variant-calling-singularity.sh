#!/bin/bash
#SBATCH -c 1
#SBATCH --constraint=skylake

bash $HOME/ecoli-analysis/run-variant-calling-singularity.sh 
