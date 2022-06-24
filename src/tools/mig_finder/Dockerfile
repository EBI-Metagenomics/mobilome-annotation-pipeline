FROM quay.io/microbiome-informatics/miniconda:4.11.0 as build

# BUILD HattCI #
ENV DEBIAN_FRONTEND noninteractive
RUN apt update && \
    apt install -y \
    python-dev \
    gcc \
    cpp \
    make \
    wget \
    zip &&\
    rm -rf /var/lib/apt/lists/*

RUN mkdir /dependencies

ENV HattCI_COMMIT="653f63278eed65bca15d3c6646d1e76ea93b2f6b"

WORKDIR /install

RUN wget https://github.com/maribuon/HattCI/archive/${HattCI_COMMIT}.zip && \
    unzip ${HattCI_COMMIT}.zip && \ 
    rm ${HattCI_COMMIT}.zip && \
    mv HattCI-"${HattCI_COMMIT}" /dependencies/hatt-ci

WORKDIR /dependencies/hatt-ci/

RUN make

# END BUILD HattCI #

RUN mamba install -c conda-forge conda-pack

RUN conda create -y -n migfinder python=2.7

ENV CONDA_DEFAULT_ENV migfinder

RUN mamba install -c bioconda -y \
    prodigal=2.6.3 \
    infernal=1.1.2 && \
    conda clean -afy

RUN . /opt/miniconda/etc/profile.d/conda.sh && conda activate migfinder && \
    pip install biopython==1.76

# Install MIG_finder #
ENV MIG_FINDER_COMMIT="b37a35dafc92254bb8eef8659efe3d95e3d86c36"

WORKDIR /install

RUN wget https://github.com/maribuon/migfinder/archive/${MIG_FINDER_COMMIT}.zip && \ 
    unzip ${MIG_FINDER_COMMIT}.zip && \
    rm ${MIG_FINDER_COMMIT}.zip && \
    mv migfinder-${MIG_FINDER_COMMIT} migfinder-src

RUN conda pack -n migfinder -o migfinder_env.tar.gz

# FINAL IMAGE #
FROM ubuntu:20.04

ENV VERSION="1.0"

LABEL maintainer="Microbiome Informatics Team www.ebi.ac.uk/metagenomics"
LABEL software="MIG_finder"
LABEL version="${VERSION}"
LABEL description="Metagenomic Integron-associated Gene finder"
LABEL website="https://github.com/maribuon/migfinder"
LABEL license="BSD-3-Clause license"

COPY --from=build /install/migfinder-src /app/migfinder
COPY --from=build /dependencies/hatt-ci /dependencies/hatt-ci
COPY --from=build /install/migfinder_env.tar.gz /app/migfinder_env.tar.gz

# HattCI dependency
RUN apt update && \
    apt install -y \
    libgomp1  &&\
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

RUN mkdir migfinder_env && \
    tar -xzf migfinder_env.tar.gz -C migfinder_env && \
    rm migfinder_env.tar.gz

SHELL ["/bin/bash", "-c"]

RUN source /app/migfinder_env/bin/activate && \
    conda-unpack

COPY run_migfinder.py /app/run_migfinder.py

ENV PATH /dependencies/hatt-ci/:/app/:/app/migfinder_env/bin:$PATH
ENV PYTHONPATH /app/migfinder:/app/migfinder_env:$PYTHONPATH

RUN ln -s /app/migfinder_env/bin/python /usr/bin/python 

# for singularity compatibility
ENV LC_ALL=C

CMD ["/bin/bash", "-c"]
