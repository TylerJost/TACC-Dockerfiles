FROM condaforge/mambaforge
RUN conda install -c conda-forge pandas
COPY runScript.py /code/runScript.py
USER root
RUN chmod +rx /code/runScript.py
ENV PATH "/code:$PATH"
