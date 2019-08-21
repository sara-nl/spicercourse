## Testing advanced features of Spider
This section includes several different features of the platform. It consists of several parts and each part contains an example that can be run independently. You can run the examples in the proposed sequence or simply pick your favourite flavor to start exploring the Spider features!

  * [High throughput data processing model](#high-throughput-data-processing-model)
  * [Software portability with containers](#software-portability-with-containers)
  * [Integration with scalable external storage](#integration-with-scalable-external-storage)
  * [Interactive analysis with Jupyter Notebooks](#interactive-analysis-with-jupyter-notebooks)
  * [Interoperability with existing platforms](#interoperability-with-existing-platforms)
* [Feedback](#4-feedback)

### High throughput data processing model

Spider is a high-throughput data-processing platform which means enabling processing of large structured data sets in short time spans. A method to achieve efficient data I/O is to split up the data processing pipelines into many parallel independent jobs where each job retrieves a chunk of data to process on the *local scratch* storage of a worker node (e.g. SSDs). This data processing model is called `embarrassingly parallel` jobs, or else the known `Grid processing model`.

Even if you never used the Grid before, this example will show you how to increase the I/O performance of your jobs by using  *local scratch* and fast network connections between the processing nodes and your input/output data storage locations. Particularly, you will:

- run a data analysis with input/output data located on your project space (on CephFS; Ceph File System)
- copy input/output data to/from a large scratch area on local SSD
- make use of the globally defined variable `$TMPDIR`

Interested? Try out the example [here](extras/tmpdir-usage.md).

### Software portability with containers

Your analysis on Spider can be run with software that was installed either by the software manager in your project or our sys admin of the system or yourself. What if you want to run the same analysis on another system(s)? Or you want to simply test some workflow on Spider but don't want to install the necessary software from scratch or mess up an existing installation? This is where containers can come in extremely handy. As you do not have admin rights on the system, we do not support Docker containers. But the good news is that we do support `Singularity` containers! In this example, you will:

- use a Singularity image with pre-installed software
- employ the container into your workflows
- run the analysis by using the software from a Singularity container

Interested? Try out the example [here](extras/singularity-usage.md).

### Integration with scalable external storage

You can download your raw data on Spider before you start your analysis. However, if you need to analyse data in excess of hundreds of TBs, your data is probably stored elsewhere, most likely on a scalable external storage system. Such a storage system for storing and retrieving huge amounts of data, is `dCache`. `dCache` consists of magnetic tape storage and hard disk storage and both are addressed by a common file system. Wouldn't it be convenient to simply download the data to be analyzed on the fly from dCache?

This can be achieved thanks to `dCache macaroons` and the high-bandwidth connection between dCache and Spider (up to 1200 Gbit/s). In this example, you will:

- make use of `dCache macaroons` within Spider
- run your analysis by fetching data on the fly directly from dCache
- evaluate options for large scale data analysis automation

Interested? Try out the example [here](extras/macaroons-usage.md).

### Interactive analysis with Jupyter Notebooks

One of the great things about Spider is that -opposite to the Grid- it allows for interactive analysis of large volumes of data. For this purpose we offer Jupyter Notebooks that can be launched on the same powerful high-throughput infrastructure of Spider. The existing Spider users can use Notebooks to run interactively some analysis with data stored/produced on Spider, or debug/prototype their pipelines with software already installed on their project spaces before submitting production runs, or even as a way to pack their work for replication and sharing with their colleagues. The Spider Notebooks is a supporting tool for all this, but not offered as a standalone service for non Spider users or training in classrooms.

Let's see how to use Jupyter Notebooks on Spider. In this example, you will:

- launch a Notebook and inspect the environment
- install packages within your Notebook or use existing software to run an analysis
- display and publish your results

Interested? Try out the example [here](extras/jupyter-usage.md).

### Interoperability with existing platforms

Recompiling your software every time you switch processing platforms or moving data around different systems is both time-consuming and makes reproducibility of your work difficult. Spider aims to be a connecting platform and our answer to the interoperability challenges between systems is: `Singularity` for software portability (see [Software portability with containers](#software-portability-with-containers)), `CVMFS/Softdrive` for software distribution and `dCache macaroons` as a user-friendly interface to large storage systems (see [Integration with scalable external storage](#integration-with-scalable-external-storage)). We will combine all this in one example that can run without modifications on HPC Cloud and Spider. In particular, we will:

- wrap our software in a Singularity container
- install the container in a central place and distribute it automagically across different systems with CVMFS/Softdrive
- get our input from dCache with macaroons (no Grid certificates)
- fetch our code from Github
- run the same analysis on multiple platforms, no vendor lock-in

Interested? Try out the example [here](extras/cloud-usage.md).

## Feedback

We hope that you enjoyed your first journey with Spider! You probably have some feedback for the platform as well as the course and we would like to hear this :)

Please fill in the survey [here]. It will take you less than 5 minutes but help us a lot!

The docuemntation for the Spider platform can be found [here](http://doc.spider.surfsara.nl/en/latest/). You can contact us with any questions by writing to helpdesk@surfsara.nl.
