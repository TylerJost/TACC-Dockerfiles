FROM continuumio/miniconda3

WORKDIR /app

# Make RUN commands use `bash --login`:
SHELL ["/bin/bash", "--login", "-c"]

# Create the environment:
COPY conda-linux-64.lock .
RUN conda create --name dev --file conda-linux-64.lock
RUN conda init bash

