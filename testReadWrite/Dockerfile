FROM condaforge/mambaforge
RUN conda install -c conda-forge pandas
USER root
COPY runScript.py /code/runScript.py
RUN chmod +rx /code/runScript.py
ENV PATH "/code:$PATH"
