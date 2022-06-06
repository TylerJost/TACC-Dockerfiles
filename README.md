# TACC-Dockerfiles

This is a repository to update and maintain dockerfiles used for single-cell sequencing analysis. 

## scanpy_development
This is a full environment for processing raw sequencing data into a log-transformed count matrix. Each portion of the `scanpy` "best-practices" paper should be possible, though some elements haven't been updated in a long time. Use of a docker environment is helpful to avoid future syntax changes. 

## sc-post
This is a post-processing development environment for single-cell sequencing through `scanpy`. It has support for most scientific computing tasks, as well as several useful machine learning packages and message passing for python. This is considered to be a template for now but will be updated to include anacoda activation support in the future.

## pod-scanpy
This is an environment that is used to run a container to develop code on through the POD. A volume is mounted which mirrors the files on Lonestar. 
 
