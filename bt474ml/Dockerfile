FROM taj159/scanpy_development:0.2
USER root
COPY runScript.py /tmp/runScript.py
RUN chmod +rx /tmp/runScript.py
ENV PATH "/tmp:$PATH"
USER $MAMBA_USER

