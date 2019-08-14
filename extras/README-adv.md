## Welcome bla bla bla

#### 1. Getting familiar with Spider - point to a welcome to spider.md and then/ or directly to our wiki

#### 2. Running a genomics pipeline

We will be using a genomics pipeline example to test some of the functionalities of Spider listed below.  The data we are going 
to use is part of a long-term evolution of E. Coli experiment led by [Richard Lenski](https://en.wikipedia.org/wiki/E._coli_long-term_evolution_experiment)
to assess adaptation in E. coli. A population was propagated for more than 50,000 
generations in a glucose-limited minimal medium. We will be working with three sample events from the Ara-3 strain of this 
experiment, one from 5,000 generations, one from 15,000 generations, and one from 50,000 generations to study how the 
population changed. Generally, the quality of raw data is assessed and data is 'trimmed'. In this example, you will download 
a small set of data that has already been trimmed and will run the variant calling workflow.

The instructions for data and software management for the examples below are already included in each example and the examples 
can be run independently. So pick your favourite feature:

1. [Analysis on local scratch space on worker nodes](https://github.com/sara-nl/spidercourse/blob/master/extras/tmpdir-usage.md)

2. [Analysis with containers](https://github.com/sara-nl/spidercourse/blob/master/extras/singularity-usage.md)

3. [Accessing data from external storage systems](https://github.com/sara-nl/spidercourse/blob/master/extras/macaroons-usage.md)

### 2. Testing advanced features with additional examples

1. Using Scientific catalogues (ND)

2. Using Jupyter notebooks (ND)

3. Infinite partition for long running jobs (RO)

4. Softdrive for software distribution (ND - if possible)

5. System performance - Pretty graphs fio example tmpdir cephfs (RO - for advisors pilot) also add massive scaling 
