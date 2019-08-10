### Using singularity containers


So far you ran your analysis with software that was installed either by the software manager or yourself. What if you want to run 
the same analysis on another system? Or you want to simply test some workflow on Spider but don't want to install the  necessary software from scratch? This is where
container can come in extremely handy. As you do not have admin rights on the system we do not support running Docker images
but we do support singularity! Let us run the same analysis but by importing software from a singularit container.

Let us first see what version of singularity is available on the system

```sh
singularity version
3.1.0-1.osgup.el7
```
In this example we will directly provide you a read-only singularity image and run the same workflow. You may first want to reset
the software environment that was set up earlier. 

```sh
cd $HOME
nano .bashrc

#Please uncomment all the lines after the following line in this file

# >>> conda initialize >>>

#Save and Exit

exit
```

Please login to Spider again and check if the previously used sofgtware is still available to you

```
fastqc -h
trimmomatic
```

This will throw errors which means thatt he software is no longer available to you. Let us now set up the scripts to use 
the Singularity containers


> **_Food for brain:_**
>
> * What does the time command do? How do you interpret the output?
> * You need to rerun the previous example with data in the project space by adding the 'time' command.
> * Does the $TMPDIR example have better performance? When is it advantageous to use it?
