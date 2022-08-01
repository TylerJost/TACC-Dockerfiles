# TACC-Dockerfiles

This is a repository to update and maintain dockerfiles used for single-cell sequencing analysis. 

## scanpy_development
This is a full environment for processing raw sequencing data into a log-transformed count matrix. Each portion of the `scanpy` "best-practices" paper should be possible, though some elements haven't been updated in a long time. Use of a docker environment is helpful to avoid future syntax changes. This is built on miniconda, so it is quite large. It has an environment called `dev` for processing data.

## sc-post
This is a post-processing development environment for single-cell sequencing through `scanpy`. It has support for most scientific computing tasks, as well as several useful machine learning packages and message passing for python. Primarily, it is less beefy than the "development" dockerfile. 

# Usage Principles
Images can currently be found at:
- https://hub.docker.com/repository/docker/taj159/scanpy_development
- https://hub.docker.com/repository/docker/taj159/sc-post
## Starting a container on the POD
The POD explicitly uses Docker. To use the images, you must start a container like so:

```
docker run -it -v /stor/work/Brock/Tyler/BT474Project:/dev/BT474Project -p 127.0.0.1:5000:8080 --name scanpy_dev taj159/scanpy_development:0.2
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







 
