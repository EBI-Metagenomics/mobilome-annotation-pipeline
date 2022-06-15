# Dockerfile for MIG_finder
FROM quay.io/microbiome-informatics/miniconda:4.11.0


ENV VERSION="1.0"

LABEL maintainer="Microbiome Informatics Team www.ebi.ac.uk/metagenomics"
LABEL software="MIG_finder"
LABEL version="${VERSION}"
LABEL description="Metagenomic Integron-associated Gene finder"
LABEL website="https://github.com/maribuon/migfinder"
LABEL license="BSD-3-Clause license"

# Dependencies
ENV DEBIAN_FRONTEND noninteractive
RUN apt update && \
    apt upgrade -y && \
    apt install -y \
        cpp \
        gcc \
        make \
        zip \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir /app /install /input /output

WORKDIR /install
RUN wget https://github.com/maribuon/HattCI/archive/653f63278eed65bca15d3c6646d1e76ea93b2f6b.zip
RUN unzip 653f63278eed65bca15d3c6646d1e76ea93b2f6b.zip \ 
    && rm 653f63278eed65bca15d3c6646d1e76ea93b2f6b.zip \
    && mv HattCI-653f63278eed65bca15d3c6646d1e76ea93b2f6b HattCI-master
WORKDIR /install/HattCI-master/
RUN make

# Install MIGfinder
WORKDIR /install
RUN wget https://github.com/maribuon/migfinder/archive/b37a35dafc92254bb8eef8659efe3d95e3d86c36.zip
RUN unzip b37a35dafc92254bb8eef8659efe3d95e3d86c36.zip \
    && rm b37a35dafc92254bb8eef8659efe3d95e3d86c36.zip \
    && mv migfinder-b37a35dafc92254bb8eef8659efe3d95e3d86c36 migfinder-master
WORKDIR /install/migfinder-master
RUN pip install -e migfinder

COPY ./run_mig.py /output

# Create the conda environment
COPY ./environment.yml /app
WORKDIR /app
RUN conda env create -f environment.yml \
   && conda clean -afy \
   && conda init bash \
   && . /root/.bashrc \
   && conda activate migfinder \
   && export PATH="/install/HattCI-master:$PATH" \
   && cd /output \
   && python ./run_mig.py 

# for singularity compatibility
ENV LC_ALL=C

