############################################################
# Dockerfile to build CRISPResso
############################################################

# Set the base image to anaconda python 3
FROM continuumio/anaconda3

# File Author / Maintainer
MAINTAINER Tony Reina/Luca Pinello 

ENV SHELL bash

RUN apt-get update && apt-get install default-jre samtools bowtie2 make gcc g++ zlib1g-dev zlib1g unzip -y 

RUN conda install biopython

RUN wget https://github.com/tonyreina/CRISPResso/archive/master.zip

RUN unzip master.zip

RUN cd CRISPResso-master && python setup.py install

ENV PATH /root/CRISPResso_dependencies/bin:$PATH

RUN rm -Rf CRISPResso-master 












