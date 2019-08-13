### Accessing data from SURfsara's dCache storage system 

You can download data on Spider before you start your analysis. However if you need to analyse data in excess of hundreds of TBs,
wouldn't it be convenient to simply download the data to be analysed on the fly? This can be achived thanks to SURFsara's large 
storage system (hard disk storage with tape backend) that provides excellent network connection (upto 1200 Gbit/s) to Spider. 
Here we will run a job where you will download the data within a job, run the analysis and push the output to your project space.

```sh
cd $HOME/ecoli-analysis
wget https://raw.githubusercontent.com/sara-nl/spidercourse/master/scripts/job-submit-variant-calling-dcache.sh

wget https://raw.githubusercontent.com/sara-nl/spidercourse/master/scripts/run-variant-calling-dcache.sh
```

Let us inspect the script that submits the job

```sh
cat job-submit-variant-calling-tmpdir.sh

#!/bin/bash
#SBATCH -c 1
#SBATCH --constraint=skylake

mkdir "$TMPDIR"/var-calling
cd "$TMPDIR"/var-calling

cp $HOME/ecoli-analysis/run-variant-calling-dcache.sh .

bash run-variant-calling-tmpdir.sh 
```

Here we first created a directory with the help of a globally defined variable $TMPDIR. This directory will be created at 
the start of the job on the local scratch space and removed when the job is done. We copy the variant calling script to this
directory and run it. Let us first run the job and while it runs we can inspect how the data transfer happens within the job.

```sh
sbatch --job-name=var-call-dcache -J 'var-call-dcache' --output=%x-%j.out job-submit-variant-calling-dcache.sh
```

Now lets inspect the variant calling script
```sh
cat run-variant-calling-tmpdir.sh

#!/bin/bash
set -e
ecolipath=$PWD

mkdir -p data/ref_genome
cp $HOME/ecoli-analysis/data/ref_genome/ecoli_rel606.fasta data/ref_genome/
ls data/ref_genome

mkdir data/trimmed_fastq_small
cp $HOME/ecoli-analysis/data/trimmed_fastq_small/*fastq data/trimmed_fastq_small/
ls data/trimmed_fastq_small

mkdir results-dcache
cd $ecolipath/results-dcache

genome=$ecolipath/data/ref_genome/ecoli_rel606.fasta

bwa index $genome

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

    bwa mem $genome $fq1 $fq2 > $sam
    samtools view -S -b $sam > $bam
    samtools sort -o $sorted_bam $bam 
    samtools index $sorted_bam
    bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam
    bcftools call --ploidy 1 -m -v -o $variants $raw_bcf 
    vcfutils.pl varFilter $variants > $final_variants
   
    done

cp -r $TMPDIR/var-calling/results-dcache $HOME/ecoli-analysis/
```

> **_Food for brain:_**
>
> * The https link from where you download the data looks rather funny. Is this a normal URL? If not, do you know what it is?
> * Is this data freely available to anyone? Try copying the link in a browser on your laptop and see what happens.
> * Dosen't look like it requires any authentication. What if your data cannot be publicly made available?

The Ecoli probably do not mind their data being public but we are all very aware of data privacy - and this is not the case only for genomic data
but research data in most domains. So how was the authnetication performed for your inout data? Try copying the link in your email
and accessing it via a browser in a week from now!

The data was shared with you with Macaroons - 
