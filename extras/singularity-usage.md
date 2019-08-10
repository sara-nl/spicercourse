### Using singularity containers


So far you ran your analysis with software that was installed either by the software manager or yourself. What if you want to run 
the same analysis on another system? Or you want to simply test some workflow on Spider but don't want to install the  necessary software from scratch? This is where
containers can come in extremely handy. As you do not have admin rights on the system, we do not support Docker containers
but we do support Singularity containers! Let us run the same analysis but by importing software from a singularity  container.

Let us first inspect what version of singularity is available on the system

```sh
singularity version
```

In this example we will directly provide you a read-only singularity image and run the same workflow. Let us get started with resetting 
the software environment that was used earlier. 

```sh
cd $HOME
nano .bashrc

#Please uncomment all the lines after the following line in this file

# >>> conda initialize >>>

#Save and Exit

exit
```

Please login to Spider again and check if the previously used software is still available to you

```
fastqc -h
trimmomatic
```

This will throw errors which means that the software is no longer available to you. Let us now download the scripts to run jobs that will use the Singularity containers

```sh
cd $HOME/ecoli-analysis
wget https://raw.githubusercontent.com/sara-nl/spidercourse/master/scripts/job-submit-variant-calling-singularity.sh

wget https://raw.githubusercontent.com/sara-nl/spidercourse/master/scripts/run-variant-calling-singularity.sh
```

Let us see how the scripts set the software environment. The job-submit-variant-calling-singularit.sh script stays the same as we do not set the software environemnt here. So, let us inspect the script that runs the analysis

```sh
cat run-variant-calling-singularity.sh

#!/bin/bash
set -e
ecolipath=$HOME/ecoli-analysis

cd $ecolipath/results

genome=$ecolipath/data/ref_genome/ecoli_rel606.fasta

singularity exec $HOME/ecoli-analysis/elixir-singularity.sif bwa index $genome

mkdir -p sam bam bcf vcf

for fq1 in $ecolipath/data/trimmed_fastq_small/*_1.trim.sub.fastq
    do
    echo "working with file $fq1"

    base=$(basename $fq1 _1.trim.sub.fastq)
    echo "base name is $base"

    fq1=$ecolipath/data/trimmed_fastq_small/${base}_1.trim.sub.fastq
    fq2=$ecolipath/data/trimmed_fastq_small/${base}_2.trim.sub.fastq
    sam=$ecolipath/results/sam/${base}.aligned.sam
    bam=$ecolipath/results/bam/${base}.aligned.bam
    sorted_bam=$ecolipath/results/bam/${base}.aligned.sorted.bam
    raw_bcf=$ecolipath/results/bcf/${base}_raw.bcf
    variants=$ecolipath/results/bcf/${base}_variants.vcf
    final_variants=$ecolipath/results/vcf/${base}_final_variants.vcf 
    
    #Align reads to reference genome
    singularity exec $HOME/ecoli-analysis/elixir-singularity.sif bwa mem $genome $fq1 $fq2 > $sam
    
    #Convert from sam to bam format
    singularity exec $HOME/ecoli-analysis/elixir-singularity.sif samtools view -S -b $sam > $bam

    #Sort the bam files    
    singularity exec $HOME/ecoli-analysis/elixir-singularity.sif samtools sort -o $sorted_bam $bam 
    
    #Calculate the read coverage of positions in the genome
    singularity exec $HOME/ecoli-analysis/elixir-singularity.sif bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam
    
    #Detect the single nucleotide polymorphisms (SNPs)
    singularity exec $HOME/ecoli-analysis/elixir-singularity.sif bcftools call --ploidy 1 -m -v -o $variants $raw_bcf 
    
    #Filter the SNPs for the final output in VCF format
    singularity exec $HOME/ecoli-analysis/elixir-singularity.sif vcfutils.pl varFilter $variants > $final_variants
   
    done
 ```

> **_Food for brain:_**
>
> * How is the software environment set up in the above script?
> * Does it matter where the singularity image is located?
> * Can you install some new software in this container? 

```
#Copy the container in your $HOME

cp /project/spidercourse/Software/elixir-singularity.sif  $HOME/ecoli-analysis

#Make sure the path to store the results in the variant calling script does not already have the results when you run the job
sbatch --job-name=var-call-singularity -J 'var-call-singularity' --output=%x-%j.out job-submit-variant-calling-singularity.sh
```

Did the analysis run successfully? You got the same results with the container without setting a single path for your software! 

Wondering how the container was built? You can find the recipe [here](https://raw.githubusercontent.com/sara-nl/spidercourse/master/extras/singularity-recipe). This container was built on Singularity version 3.2.1-1 (on macOS Mojave running Vagrant). Using containers is simple but you should bear in mind that it should be properly built with all the required dependencies, and should be updated regularly and tested in the runtime environment. 
