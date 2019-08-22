## Interoperability with existing platforms 

In this section you will run the same E. coli analysis example to test the interoperability features of Spider but on another SURFsara platform- the HPC Cloud.

In this example we have pre-configured a Virtual Machine with the necessary components (Data and Software) to integrate the analysis with Spider. In order to access the VM, ask the instructors for your credentials. 

### Run the analysis

* Login to the pre-configured VM on HPC Cloud (ask for your username/password if you still don't have this). Once logged in download the analysis script and execute it:

```sh
ssh username@spider.usersupport-cloud.surf-hosted.nl 
cd $HOME
mkdir ecoli-analysis-cloud
cd ecoli-analysis-cloud/
wget https://raw.githubusercontent.com/sara-nl/spidercourse/master/scripts/run-variant-calling-cloud-adv.sh
chmod u+x run-variant-calling-cloud-adv.sh
./run-variant-calling-cloud-adv.sh
```

* While your analysis is running inspect the script `run-variant-calling-cloud-adv.sh`. Where is your input data fetched from? Where is your software installed? How is the Software ported to this VM on HPC Cloud?

> Food for brain:
> - Run the same analysis on Spider. How would you submit the `run-variant` script to the cluster? Do you have to make any changes to your software or input data paths?  
> - Think of other platforms to port the same analysis, e.g. Cartesius. Sketch your solution.

Curious how this was set up? Wanna give it a go on your own laptop?
 - Install and configure CVMFS for software distribution- follow the instructions [here](http://doc.grid.surfsara.nl/en/latest/Pages/Advanced/softdrive_on_laptop.html#softdrive-on-laptop)
 - Install Singularity - follow the instructions [here](https://sylabs.io/guides/3.0/user-guide/installation.html) and install from source Singularity v3.1.0
