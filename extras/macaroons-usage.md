### Accessing data from SURfsara's dCache storage system 

You can download your raw data on Spider before you start your analysis. However, if you need to analyse data in excess of hundreds of TBs,
wouldn't it be convenient to simply download the data to be analysed on the fly? This can be achieved thanks to SURFsara's large 
storage system (hard disk storage with tape backend) that provides excellent network connection (upto 1200 Gbit/s) to Spider. 
Here we will run a job where you will download the data directly on the worker node within the job, run the analysis and push the output to your home directory.

```sh
cd $HOME/ecoli-analysis
wget https://raw.githubusercontent.com/sara-nl/spidercourse/master/scripts/job-submit-variant-calling-dcache.sh

wget https://raw.githubusercontent.com/sara-nl/spidercourse/master/scripts/run-variant-calling-dcache.sh
```

Let us inspect the script that submits the job

```sh
cat job-submit-variant-calling-dcache.sh

#!/bin/bash
#SBATCH -c 1
#SBATCH --constraint=skylake

mkdir "$TMPDIR"/var-calling
cd "$TMPDIR"/var-calling

cp $HOME/ecoli-analysis/run-variant-calling-dcache.sh .

bash run-variant-calling-dcache.sh 
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

cd data/
curl https://webdav.grid.surfsara.nl:2880/?authz=MDAxY2xvY2F0aW9uIE9wdGlvbmFsLmVtcHR5CjAwMThpZGVudGlmaWVyIDVMdFI5S29QCjAwMzJjaWQgaWQ6NDM2MzI7NDEzODUsNDQ0MzYsNDI1MjksMzAwMTM7bWFpdGhpbGsKMDAyOGNpZCBiZWZvcmU6MjAxOS0wOS0xMlQxMDoxMzoyNy42NzVaCjAwNWFjaWQgcm9vdDovcG5mcy9ncmlkLnNhcmEubmwvZGF0YS9sc2dyaWQvU1VSRnNhcmEvc3BpZGVyY291cnNlL3RyaW1tZWRfZmFzdHFfc21hbGwudGFyCjAwMWZjaWQgYWN0aXZpdHk6RE9XTkxPQUQsTElTVAowMDJmc2lnbmF0dXJlIGL5MfchTf7sH1Ela025OBtIiYmsB3LAbutPyTbgW73yCg --output trimmed_fastq_small.tar
tar xvf trimmed_fastq_small.tar

mkdir $ecolipath/results-dcache
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
but research data in most domains. So how was the authnetication performed for your input data? 

The data was shared with you with Macaroons - these are bearer tokens that you can use to authorize someone to dwonload/upload/delete data stored on dCache. These macaroons can be used with clients that can support bearer tokens (e.g., curl, Rclone). For this exercise a macaroon was created with certain restrictions (called as caveats) on the lifetime of the macaroon, the IP address you can use the macaroon from, the file that you can access, etc. Depending on who you want to share the data with, for how long and from which systems, these caveats can be adjusted. 

You can check the status of the job and inspect the output log file (even if the job is not completed).

```sh
squeue -u $USER
cat var-call-dcache-jobid.out #replace the jobid with your jobid 
```

You downloaded the input data on the fly (if you have good netowrk connectivity to the storage system from Spider) within each job without the hassle of downloading all the data to Spider. This is particularly handy if you want to automate large scale data analysis. In this example we saved the results locally on Spider, but you can also push the output to dCache or another external storage system. dCache also supports username/password authentication and certificate based proxy authentication.   
