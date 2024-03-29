# TACC-Dockerfiles

This is a repository to update and maintain dockerfiles used for single-cell sequencing analysis. 

If you are curious about what Docker is, check out the always helpful Fireship video: https://www.youtube.com/watch?v=Gjnup-PuquQ

## scanpy_development
This is a full environment for processing raw sequencing data into a log-transformed count matrix. Each portion of the `scanpy` "best-practices" paper should be possible, though some elements haven't been updated in a long time. Use of a docker environment is helpful to avoid future syntax changes. This is built on miniconda, so it is quite large. It has an environment called `dev` for processing data.

## sc-post
Don't use this one, it is under quarantine. 

# Usage Principles
Images can currently be found at:
- https://hub.docker.com/repository/docker/taj159/scanpy_development

To download with Docker:
```
docker pull taj159/scanpy_development:0.5
```
- https://hub.docker.com/repository/docker/taj159/sc-post

To download with Docker:
```
docker pull taj159/sc-post:0.6
```

## Starting a container on the POD
The POD explicitly uses Docker. To use the images, you must start a container like so:

```
docker run -it -v /stor/work/Brock/Tyler/BT474Project:/dev/BT474Project -p 127.0.0.1:5000:8888 --name scanpy_dev taj159/scanpy_development:0.5
```

Breaking this down, we are running this image with the following terms:
| Flag      | Description |
| ----------- | ----------- |
| -it       | Runs interactively       |
| -v   | Mounts a local volume inside the container<br />(`local directory:remote directory`)        |
| -p        | Exposes a port (`localhost:localport:remote port`) |
| --name | The name of the container|

Finally, we end with the name of the image. We can start the container after exiting using:

```
docker start -i scanpy
```

### Attaching to Containers in VS Code
Once you have a running container, you can choose to start a Jupyter notebook if you wish. If you prefer to work in an editor such as VS Code, you must attach to a *running* container. To do so, first install the [Docker extension](https://code.visualstudio.com/docs/containers/overview). You can then start the extension, select the container, then select "Attach Visual Studio Code". 

![](./media/containerAttach.png).  

*Note*: You will have to install extensions for VS Code inside the container that you might already have, more instructions can be found [here](containerAttach.png). 
## Starting a container through TACC
TACC uses singularity instead of Docker, so things are a bit different. To use singularity, first call the module with `module load tacc-singularity`. 

To begin, first download the image using singularity. In this case, we would call `singularity pull docker://taj159/scanpy_development:0.2`. This will download a `.sif` file that we will use for the remainder. 

Singularity is actually quite nice because it (currently) will expose all directories and ports. For our image, we can call:

```
singularity run ./scanpy_development_0.2.sif
```
From here, we can proceed as normal. If we were running this from a slurm file, we would need to give it extra commands, such as:

```
singularity exec ./scanpy_development_0.2.sif conda run --no-capture-output -n dev python <script>
```

This would activate the conda environment `dev` on start and run our desired script.

# TACC/Docker/Singularity Cheatsheet
## Docker/Singularity Commands
```
docker image ls #List image
docker ps –a #List all containers
docker build -t <username>/<image name>:<tag> .
docker image push <username>/<image name>:<tag>
docker run -it -v <folder> -p <port forward> --name <container> <image> 

singularity pull docker://<username>/<image>
singularity exec <file>.sif
singularity shell <file>.sif
```
## TACC Commands
```
module list #Lists loaded TACC modules
module load <module name> #Load your module
sbatch submission.slurm #Submit your job
sacct -j <myid> --format=Elapsed #Time to complete job
idev –m 40 #Starts interactive session
squeue –u <username> #Check job status
module spider <term> # Searches for term
idev –m <Num Minutes> -p <queue name>
```
## TACC Slurm Template
```
#!/bin/bash
#SBATCH -J test_job                   		# Job name
#SBATCH -o test_job.%j.out        		# Name of stdout output file (%j expands to jobId)
#SBATCH -p normal   	      			# Queue name (you use the development queue to test jobs that take less than 2 hours)
#SBATCH -N 1                  			# Total number of nodes requested
#SBATCH -n 1                 			# Total number of threas tasks requested (128 per node)
#SBATCH -t 01:30:00           			# Run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=email@utexas.edu   	        # Address email notifications
#SBATCH --mail-type=all				# Email at begin and end of job
#SBATCH -A DMS21043      			# <-- Allocation name to charge job against

module load tacc-singularity
singularity run docker://<username>/<image name>:<tag> #Run default command
singularity exec docker://<username>/<image name>:<tag> conda run --no-capture-output -n <environment> python <script> #Run arbitrary script using the given environment
```
# Other Notes for Effective Use
## Singularity/Docker/Anaconda Peculiarities
Space is a crucial portion of any docker file. The full scanpy development docker file rests at about 5 GB (for reference the maximum that you can upload is around 10GB after compression). This is fairly hefty for Docker, but can be reduced by using [micromamba](https://github.com/mamba-org/micromamba-docker). However, this might not be ideal for you.  

If you used micromamba, you are free to use `singularity run <.sif file> <commands>`. However, if you gain shell access either through `singularity shell` or `singularity run`, you won't be able to access the anaconda environment. To remedy this, input:
```
eval "$(micromamba shell hook --shell=bash)" && micromamba activate
``` 
On startup and your environment should activate. 

If you are running an anaconda on TACC and have a miniconda-created conda environment in your Docker file, you may have to specify which conda you want to run. This tends to get messy quickly. 
