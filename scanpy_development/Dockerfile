FROM continuumio/miniconda3

WORKDIR /app

# Make RUN commands use `bash --login`:
SHELL ["/bin/bash", "--login", "-c"]

# Create the environment:
COPY environment.yml .
RUN conda env create -f environment.yml

RUN conda init bash

RUN conda env export --name dev > env.lock
