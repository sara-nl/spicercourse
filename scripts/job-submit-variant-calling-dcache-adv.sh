#!/bin/bash
#SBATCH -c 1
#SBATCH --constraint=skylake

mkdir "$TMPDIR"/var-calling
cd "$TMPDIR"/var-calling

cp $HOME/ecoli-analysis-dcache/run-variant-calling-dcache-adv.sh .

export PATH="/project/spidercourse/Software/ecoli-analysis-software/miniconda2/bin:$PATH"

bash run-variant-calling-dcache-adv.sh 
