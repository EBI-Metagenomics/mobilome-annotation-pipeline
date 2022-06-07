# Dockerfile for ICE_finder
FROM quay.io/microbiome-informatics/miniconda:4.11.0


ENV VERSION="1.0"

LABEL maintainer="Microbiome Informatics Team www.ebi.ac.uk/metagenomics"
LABEL software="ICE_finder"
LABEL version="${VERSION}"
LABEL description="A tool for the detection of ICEs/IMEs of bacterial genomes"
LABEL website="https://bioinfo-mml.sjtu.edu.cn/ICEfinder/ICEfinder.html"
LABEL license="The Perl-based local version of ICEfinder is available by contact Dr Hong-Yu (hyou@sjtu.edu.cn)"

# Dependencies
ENV DEBIAN_FRONTEND noninteractive
RUN apt update && \
    apt upgrade -y && \
    apt install -y \
        cpp \
        curl \
        gcc \
        git \
        make \
    && rm -rf /var/lib/apt/lists/*

RUN conda config --add channels defaults && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda

RUN mamba install -y \
    blast=2.12.0 \
    emboss=6.6.0 \
    perl-bioperl \
    && conda clean -afy

RUN pip install gdown

RUN cpanm Linux::Distribution Getopt::Std

# Set environment
RUN export PERL5LIB="/opt/miniconda/lib/perl5/site_perl/5.22.0/Bio"

## Getting ICE_finder
# The local version of ICEfinder is located in a googleDrive folder
# This is a temporary solution to get the tool available
# Downloading the tar file
RUN mkdir /install /input /output /work
WORKDIR /install
RUN gdown --id 1RlkjErOg08W-MUrO7cLUm_G3FLB2Yq7g
RUN tar -xvzf ICEfinder_linux.tar.gz ICEfinder_linux 
RUN rm ICEfinder_linux.tar.gz

# Fixing tools location
WORKDIR /install/ICEfinder_linux/tools
RUN rm seqret transeq blastn blastp
RUN ln -s /opt/miniconda/bin/seqret seqret \
    && ln -s /opt/miniconda/bin/transeq transeq \
    && ln -s /opt/miniconda/bin/blastn blastn \
    && ln -s /opt/miniconda/bin/blastn blastp

# Copy data to analyse
COPY ./gbk_files/input.list /install/ICEfinder_linux/
COPY ./gbk_files/*.gbk /install/ICEfinder_linux/gbk/

WORKDIR /install/ICEfinder_linux/

# for singularity compatibility
ENV LC_ALL=C

