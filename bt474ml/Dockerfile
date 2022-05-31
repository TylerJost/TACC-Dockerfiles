FROM taj159/sc-post:0.1
# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "sc-post", "/bin/bash", "-c"]

# Demonstrate the environment is activated:
RUN echo "Make sure scanpy is installed:"
RUN python -c "import scanpy"
# The code to run when container is started:
COPY runScript.py /code/runScript.py
RUN chmod +rx /code/runScript.py
ENV PATH "/code:$PATH"
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "sc-post", "python", "/code/runScript.py"]
