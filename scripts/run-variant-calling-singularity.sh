#!/bin/bash
set -e
ecolipath=$HOME/ecoli-analysis

cd $ecolipath
mkdir $SLURM_JOBID-results
cd $SLURM_JOBID-results

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
    sam=$ecolipath/$SLURM_JOBID-results/sam/${base}.aligned.sam
    bam=$ecolipath/$SLURM_JOBID-results/bam/${base}.aligned.bam
    sorted_bam=$ecolipath/$SLURM_JOBID-results/bam/${base}.aligned.sorted.bam
    raw_bcf=$ecolipath/$SLURM_JOBID-results/bcf/${base}_raw.bcf
    variants=$ecolipath/$SLURM_JOBID-results/bcf/${base}_variants.vcf
    final_variants=$ecolipath/$SLURM_JOBID-results/vcf/${base}_final_variants.vcf 
    
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
