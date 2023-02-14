############################################################
# Dockerfile to build CRISPResso
############################################################

FROM mambaorg/micromamba
# File Author / Maintainer
MAINTAINER Tony Reina/Luca Pinello

WORKDIR /app

COPY *.yml /app
COPY *.py /app
COPY *.txt /app
COPY LICENSE /app
COPY README.md /app
COPY CRISPResso /app/CRISPResso

SHELL ["/bin/bash", "-c"]

RUN micromamba install -y -n base -f environment.yml && micromamba clean --all --yes

SHELL ["mamba", "run", "--name", "base", "/bin/bash", "-c"]
