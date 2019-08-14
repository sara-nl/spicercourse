#!/bin/bash
#SBATCH -c 1
#SBATCH --constraint=skylake

mkdir "$TMPDIR"/var-calling
cd "$TMPDIR"/var-calling

cp $HOME/ecoli-analysis-tmpdir/run-variant-calling-tmpdir-adv.sh .

export PATH="/project/spidercourse/Software/ecoli-analysis-software/miniconda2/bin:$PATH"

time bash run-variant-calling-tmpdir.sh 
