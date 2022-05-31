# Dockerfile.just-mamba
FROM condaforge/mambaforge
COPY environment.yml .
RUN mamba env create -f environment.yml
