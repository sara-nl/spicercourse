# Spider project spaces and roles

1. [Spider project environment](#spider-spaces)
2. [Data management](#spider-dm)
3. [Software management](#spider-sm)
4. [Data inspection and cleanup](#data-cleanup)

### <a name="spider-spaces"></a> 1. Spider project environment

#### 1.1 Project environment

Each project on Spider gets its own project space. Let us get familiar with this environment :

 ```sh
 ls /project 
 ```
> **_Food for brain:_**
>
> * Do you know what project you belong to? Hint: Combine the information from the commands "whoami" and "ls /project"
> * What are each of the project directories for? What is public/private? Do you have read and write permissions on all spaces? Hint: Does the command "echo "Hello world" > /project/yourproject/$USER-hello-world.txt" create the file for you? What about the directories in /project/yourproject/?

### <a name="spider-dm"></a> 2. Data management

#### 2.1 Data manager Role :information_desk_person:

 ```sh
 id $USER
 ```

> **_Food for brain:_**
>
> * What is the role of a data manager? Do you have that role? If not, do you know who the data manager is?

Let us use some handy commands that can tell you which 'groups' different users belong to.

 ```sh
 getent group spidercourse-data
 getent group spidercourse-user
 ```
 Which group does your username belong to? And which one does it not belong to?
 
 ```sh
 mkdir /project/spidercourse/Data/mydata
 ```
 
What happened? As you are not a data manager :information_desk_person: (so not a part of spider_data group) but a user (part of spider_user group), you do not have write permissions in the project's Data directory. However, all users in this project have read permissions. This means data is managed by a 'data manager' and can be shared within the project without worrying about regular users accidentally deleting/overwriting the data for the project.

#### 2.2 Data download

Let us download some data that we will use later to run jobs on Spider. The data we are going to use is part of a long-term evolution experiment led by [Richard Lenski](https://en.wikipedia.org/wiki/E._coli_long-term_evolution_experiment) to assess adaptation in E. coli. A population was propagated for more than 50,000 generations in a glucose-limited minimal medium. We will be working with three sample events from the Ara-3 strain of this experiment, one from 5,000 generations, one from 15,000 generations, and one from 50,000 generations to study how the population changed. 

Let us download the paired-end data from [European Nucleotide Archive](https://www.ebi.ac.uk/ena). The specific instructions for Data manager or regular users are marked accordingly, the generic instructions should be followed by everyone. 

:information_desk_person: As Data manager

 ```sh
 mkdir -p /project/spidercourse/Data/ecoli-analysis/data/untrimmed_fastq/
 cd /project/spidercourse/Data/ecoli-analysis/data/untrimmed_fastq/
 ```
 
:construction_worker: As regular user
 
 ```sh
 mkdir -p $HOME/ecoli-analysis/data/untrimmed_fastq/
 cd $HOME/ecoli-analysis/data/untrimmed_fastq/
 ```

Let us now download the data

 ```sh
 curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_1.fastq.gz
 curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_2.fastq.gz
 curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_1.fastq.gz
 curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_2.fastq.gz
 curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_1.fastq.gz
 curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_2.fastq.gz 
 
 #Sometimes the ftp server is rather slow. In this case to save tim e you can copy the data from a local folder where we have downloaded the files for you already
 cp /project/spidercourse/Data/ecoli-analysis/data/untrimmed_fastq/*fastq.gz .
 
 ```

Let us also download the reference genome for E. coli REL606.

:information_desk_person: As Data manager

 ```sh
 cd /project/spidercourse/Data/ecoli-analysis/data
 mkdir ref_genome
 cd ref_genome
 ```
 
:construction_worker: As regular user
 
 ```sh
 
 cd $HOME/ecoli-analysis/data
 mkdir ref_genome
 cd ref_genome
 ```
 
 ```sh
 curl -L -o ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
 gunzip ecoli_rel606.fasta.gz
 ```

### <a name="spider-sm"></a> 3. Software management

#### 3.1 Software manager Role :woman:

 ```sh
 id $USER
 ```
 
> **_Food for brain:_**
>
> * Are you a software manager? If not, do you know who is the software manager? 

 ```sh
 getent group spidercourse-sw
 getent group spidercourse-user
 mkdir /project/spidercourse/Software/mysoftware
 ```
 
What happened? As you are not a Software manager :woman: (so not a part of the spider_sw group), you do not have write permissions in the project's Software directory. However, all users still have read permissions. This means that you or one of your colleagues in this role can install complicated software and the dependencies, maintain it, and all project users can use that uniformly. This makes life easier and is crucial for reproducability of your results. Also, apart from project wide software, individual users can install their own software (or different versions) as you shall see below. 

#### 3.2 Miniconda installation

We will use Miniconda; it is a package manager. We shall first install miniconda2 and then proceed to the installation of individual tools. The specific instructions for Software manager or regular users are marked accordingly, the generic instructions should be followed by everyone.  

:woman: As Software manager
 
 ```sh
 cd /project/spidercourse/Software/ 
 ```
 
 :construction_worker: As regular user
 ```sh
 cd $HOME
 ```
 
 ```sh
 wget https://repo.continuum.io/miniconda/Miniconda2-4.6.14-Linux-x86_64.sh
 bash Miniconda2-4.6.14-Linux-x86_64.sh
 ```

It will ask you 

 ```sh
 Do you accept the license terms? [yes|no]
 yes # please answer yes to proceed further
 ```

It will ask you for an installation path. 

:woman: As Software manager
 
 ```sh
  /project/spidercourse/Software/ecoli-analysis-software/miniconda2 
 ```
 
 :construction_worker: As regular user
 ```sh
 $HOME/ecoli-analysis-software/miniconda2 
 ```
 
 It will then ask you 
 
 ```sh
 Do you wish the installer to initialize Miniconda2 by running conda init? [yes|no]
 yes # please answer yes
 ...
 ...
 ...
 #If you'd prefer that conda's base environment not be activated on startup, set the auto_activate_base parameter to false: 
 #conda config --set auto_activate_base false 
 
 #Thank you for installing Miniconda2!
 
 #This shows that you successfully installed miniconda2. For these changes to take effect, please logout and login again.
 exit 
 ```

We will see later how Software installed in the Software project space can be used by all users. If you installed it in your $HOME directory, other users in the project cannot access it. 

Login again to Spider and inspect what environment variables have been set up

 ```sh
 cat $HOME/.bashrc
 ```

#### 3.3 Variant calling tools installation

Follow the further instructions for the installation of individual tools

:woman: As Software manager

 ```sh
 cd /project/spidercourse/Software/ecoli-analysis-software 
 ```
 
:construction_worker: As regular user

 ```sh
 cd $HOME
 ```
 
 ```sh
 conda install -c bioconda fastqc=0.11.7=5

 conda install -c bioconda trimmomatic=0.38=0

 conda install -c bioconda bwa=0.7.17=ha92aebf_3

 conda install -c bioconda samtools=1.9=h8ee4bcc_1

 conda install -c bioconda bcftools=1.8=h4da6232_3 
 ```
 
> **_Food for brain:_**
>
> * Do you know where the above pacakages are installed? How do you test if the installation was successful? #Hint - try running one of the above installed software

```sh
fastqc -h       (fastqc runs quality control checks on raw sequence data and with the -h flag it displays the help for this task)

trimmomatic     (data trimming tool)

bwa             (maps DNA sequences against reference genome)

samtools        (utility for manipulating data in the SAM format)

bcftools        (variant calling tool)
```

For the last two tools you can see that a library is missing! This will very often be the case that not every default installation will run on every system - some dependencies might be missing and would require troubleshooting. If the Software manager resolves these issues, users can freely use the software and avoid hassles of resolving software dependencies. In this case, the Softare manager has already resolved this for you and instead of troubleshooting you can already start using this installation by performing the following steps:

```sh
nano $HOME/.bashrc

#In the conda initialize set up replace all the paths 

# >>> conda initialize >>>

replace the following path in your file 

$HOME/ecoli-analysis-software/miniconda2/bin/conda

with the following path instead

/project/spidercourse/Software/ecoli-analysis-software/miniconda2/bin/conda

#Similarly change the other respective paths directing to the Software in project space

save the changes

exit 
```
Now you are using the software environment that was set up by the Software manager. Login to Spider again.

```ss
echo $PATH

samtools
bcftools
```
As you can see the error is resolved and you can proceed to running the analysis. Just so you know in this particular case it was a simple solution, the missing library was downloaded to Spider and copied as follows "cp libcrypto.so.1.0.0 ecoli-analysis-software/miniconda2/lib/"

### <a name="data-cleanup"></a> 4. Data inspection and cleanup

Now that you have the raw data and the software installed, we will assess the quality of the sequence reads contained in our fastq files and run filtering. Pay attention to the changes you need to make to the scripts.

:information_desk_person: As data manager

```sh
cd /project/spidercourse/Data/ecoli-analysis/
```

:construction_worker: As a regular user

```sh 
cd $HOME/ecoli-analysis
```

```
wget https://raw.githubusercontent.com/sara-nl/spidercourse/master/scripts/job-submit-datatrimming.sh
 
cat job-submit-datatrimming.sh

#!/bin/bash
#SBATCH -c 1
#SBATCH --constraint=skylake

#As data manager uncomment the following line (remove #)
#bash /project/spidercourse/Data/ecoli-analysis/data_qc.sh 

#As a regular user uncomment the following line (remove #)
#bash $HOME/ecoli-analysis/data_qc.sh 
```

Let us download the data_qc.sh script inspect what steps we follow in the data trimming:

```sh
wget https://raw.githubusercontent.com/sara-nl/spidercourse/master/scripts/data_qc.sh

cat data_qc.sh 

#!/bin/bash
set -e

#if you are a data manager, uncomment the following line (remove the #)
#ecolipath=/project/spidercourse/Data/ecoli-analysis

#if you are a regular user, uncomment the following line (remove the #)
#ecolipath=$HOME/ecoli-analysis

mkdir -p $ecolipath/data/fastqc_untrimmed_reads

cd $ecolipath/data/fastqc_untrimmed_reads/

echo "Running FastQC ..."
fastqc $ecolipath/data/untrimmed_fastq/*.fastq* -o ./ 

cd $ecolipath/data/untrimmed_fastq
cp /project/spidercourse/Software/ecoli-analysis-software/miniconda2/pkgs/trimmomatic-0.38-0/share/trimmomatic-0.38-0/adapters/NexteraPE-PE.fa .
echo "Running trimmomatic"
for infile in *_1.fastq.gz
do
   base=$(basename ${infile} _1.fastq.gz)
   trimmomatic PE ${infile} ${base}_2.fastq.gz \
               ${base}_1.trim.fastq.gz ${base}_1un.trim.fastq.gz \
               ${base}_2.trim.fastq.gz ${base}_2un.trim.fastq.gz \
               SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 
done

mkdir $ecolipath/data/trimmed_fastq	
mv *.trim* $ecolipath/data/trimmed_fastq

echo "done"
```
Now let's submit the job

```sh
sbatch --job-name=data-trim -J 'data-trim' --output=%x-%j.out job-submit-datatrimming.sh
squeue -u $USER
```
This job will create trimmed data file and files in .html format that can be used to display the summary report.

> **_Food for brain:_**
>
> * The sbatch command you used here is different than for the previous job. What do the additional flags do? Hine: Check your jobname with "squeue -u $USER" and the stdout filename (e.g., earlier it was slurm-jobid.out)

In this job, we have combined the creation of files to assess quality control as well as trim the data. Typically you would first assess the data and then perform the trimming. The files generated by fastqc tool create .html files that can be used to visualize for assessing the data. These files have already been created and copied at the following path that you can inspect while your job runs:

```sh
cd /project/spidercourse/Public/
ls ecoli-analysis
```

The Public directory is yet another project space that can be accessed by not only all Spider users but is open to anyone.  This public space is exposed to the web under the domain: https://public.spider.surfsara.nl/project/[PROJECTNAME]/. To view the files go to https://public.spider.surfsara.nl/project/spidercourse/ link from the browser on your laptop. You can view one of the files in the ecoli-analhysis folder and [study](https://datacarpentry.org/wrangling-genomics/02-quality-control/index.html) the quality of the data we have.

The Public project space makes data sharing very easy and fast, particularly with your colleagues who are not a part of the project on Spider. 

#### Acknowledgements 
This example was adopted from https://datacarpentry.org/wrangling-genomics/ 

 
