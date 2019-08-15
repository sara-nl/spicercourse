## Testing advanced features of Spider

### 1. Running E-coli analysis with advanced features

If you already ran the earlier example with [data and software management](https://github.com/sara-nl/spidercourse/blob/master/demo-spider-roles.md) and [data analysis](https://github.com/sara-nl/spidercourse/blob/master/run-spider-jobs.md), you are familiar with project spaces and 
various roles in a project and may start running the following examples directly. 

If you skipped these, here is a quick explanation of what analysis you will run in this section:
- The data in use is part of a long-term evolution experiment 
led by Richard Lenski to assess adaptation in E. coli. A population was propagated for more than 50,000 generations in a 
glucose-limited minimal medium. 
- We will be working with three sample events from the Ara-3 strain of this experiment, 
one from 5,000 generations, one from 15,000 generations, and one from 50,000 generations to study how the population changed.
- You will download a small set of data that has already been trimmed and will run the variant calling workflow.

The instructions for data and software management for the examples below are included in each example and the examples can be run independently.


1. [Analysis on local scratch space on worker nodes](https://github.com/sara-nl/spidercourse/blob/master/extras/tmpdir-usage.md)

2. [Analysis with containers](https://github.com/sara-nl/spidercourse/blob/master/extras/singularity-usage.md)

3. [Accessing data from external storage systems](https://github.com/sara-nl/spidercourse/blob/master/extras/macaroons-usage.md)

### 2. Testing advanced features with additional examples

1. Using Scientific catalogues (ND)

2. Using Jupyter notebooks (ND)

3. Infinite partition for long running jobs (RO)

4. Softdrive for software distribution (ND - if possible)

5. creare SC (RO+MK)

5. System performance - Pretty graphs fio example tmpdir cephfs (RO - for advisors pilot) also add massive scaling 

6. analysis on HPC Cloud VM (ND)
