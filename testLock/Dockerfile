# Dockerfile

# Build container
FROM continuumio/miniconda:latest as conda

ADD conda-linux-64.lock /locks/conda-linux-64.lock
RUN conda create -p /opt/env --copy --file /locks/conda-linux-64.lock

# Primary container

FROM gcr.io/distroless/base-debian10

COPY --from=conda /opt/env /opt/env
