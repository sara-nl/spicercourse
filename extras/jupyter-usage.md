## Interactive analysis with Jupyter Notebooks

Let's start our first Spider notebook! To make it easy, we have created a tool to allow users 
launch Jupyter Notebooks. This is called `startnotebook` and can be invoked from the Spider User 
Interface (UI).

* Login to Spider and inspect the `startnotebook`  tool:

```sh
$ startnotebook -h
```

As you see `startnotebook` accepts the following parameters:

| option | description |
| ------ | ----------- |
|`--flavor FLAVOR` | The FLAVOR determines the amount of resources assigned to the Jupyter notebook. Currently the  The available flavours are `small` (1core/8GBmem) and `medium` (4cores/32GBmem). If not specified, the deafult is `small`. |
|`--check` | This is a dry run to inspect the resources to be used by the notebook. It doesn't launch a notebook. |
|` --kernel-timeout KERNEL_TIMEOUT` | Idle time (s) before the kernel is killed. This counts when a kernel is idle. If not specified, the deafult is 10 min. |
|`--notebook-timeout NOTEBOOK_TIMEOUT`| Time (s) before the notebook is killed. This counts when no kernels are running. If not specified, the deafult is 10 min. |
|`--name JOBNAME`| The JOBNAME is the name you want to give to your Notebook folder. If not specified, a random name is created as 'notebook-[random string]' |
|` --verbose`| Verbose mode to display output from background processes that are triggered. | 

* Start a notebook with name 'my-awesome-research`

```sh
$ startnotebook --name my-awesome-notebook
# Starting Notebook.
# Determining Notebook url......
# Notebook sucessfully launched!
# Access your notebook using the following parameters:

	url : http://145.38.255.23:8827
	pass: OdL@@F>yQ:nQ
```

> What happened? Behind the scenes, `startnotebook` submits a job to the SLURM cluster. The job executes a Singularity container, containing the Jupyter notebook. Once the job starts running the tool retrieves the host and port where the Jupyter notebook is running and provides you with the corresponding URL, and credentials to allow you access the notebook from your laptop.

* Open a browser in your laptop and copy the provided URL. Hit 'Enter' and paste the provided password. Did you successfully log in? Can you start a new 'Python 3' notebook?

* Check hidden configurtion files created locally in your Spider home account. You should see something like:

```sh
ls -la my-awesome-notebook/
# drwxrwxr-x .config
# drwxr-xr-x .ipynb_checkpoints
# drwxr-xr-x .ipython
# drwxrwxr-x .jupyter
# drwxr-xr-x .local
# -rw-r--r-- 1 Untitled.ipynb
```



