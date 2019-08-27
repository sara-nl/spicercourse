## Interoperability with existing platforms 

In this section you will run the same E. coli analysis example to test the interoperability features of Spider but on another SURFsara platform- the HPC Cloud.

In this example we have pre-configured a Virtual Machine with the necessary components (Data and Software) to integrate the analysis with Spider. In order to access the VM, ask the instructors for your credentials. 

### Run the analysis

* Login to the pre-configured VM on HPC Cloud (ask for your username/password if you still don't have this). Once logged in download the analysis script and execute it:

```sh
ssh username@spider.usersupport-cloud.surf-hosted.nl 
mkdir -p /cloud-analysis/$USER/ecoli-analysis-cloud
cd /cloud-analysis/$USER/ecoli-analysis-cloud
wget https://raw.githubusercontent.com/sara-nl/spidercourse/master/scripts/run-variant-calling-cloud.sh
chmod u+x run-variant-calling-cloud.sh
./run-variant-calling-cloud.sh
```

* While your analysis is running inspect the script `run-variant-calling-cloud.sh`.  

You just ran the same pipeline on another system!

> Food for brain:
> - Where is your input data fetched from?
> - Where is your analysis running?
> - Where is your software installed? 
> - How is the Software ported to this VM on HPC Cloud?
> - Think of other platforms to port the same analysis, e.g. Cartesius. Sketch your solution.

Curious how this was set up? Wanna give it a go on your own laptop?
 - Install and configure CVMFS for software distribution- follow the instructions [here](http://doc.grid.surfsara.nl/en/latest/Pages/Advanced/softdrive_on_laptop.html#softdrive-on-laptop)
 - Install Singularity - follow the instructions [here](https://sylabs.io/guides/3.0/user-guide/installation.html) and install from source Singularity v3.1.0
