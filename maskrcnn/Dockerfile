FROM tensorflow/tensorflow:2.0.3-gpu

MAINTAINER Tyler Jost <tyler_jost@utexas.edu>

RUN apt-get update && \
    apt-get install -y curl build-essential libpng12-dev libffi-dev  && \
    apt-get clean && \
    rm -rf /var/tmp /tmp /var/lib/apt/lists/*

RUN curl -sSL -o installer.sh https://repo.anaconda.com/archive/Anaconda3-2018.12-Linux-x86_64.sh && \
    bash /installer.sh -b -f && \
    rm /installer.sh

ENV PATH "$PATH:/root/anaconda3/bin"

