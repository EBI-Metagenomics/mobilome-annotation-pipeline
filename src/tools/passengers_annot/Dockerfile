# Dockerfile for blast annotation of integron passenger genes
FROM quay.io/microbiome-informatics/miniconda:4.11.0


ENV VERSION="2.12.0"

LABEL maintainer="Microbiome Informatics Team www.ebi.ac.uk/metagenomics"
LABEL software="Blast"
LABEL version="${VERSION}"
LABEL description="Integron passenger genes blast annotation"
LABEL webpage="https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download"
LABEL license="https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/scripts/projects/blast/LICENSE"

# Dependencies
ENV DEBIAN_FRONTEND noninteractive
RUN apt update && \
    apt upgrade -y && \
    apt install -y \
        cpp \
        dos2unix \
        gcc \
        git \
        gzip \
        make \
        zip \
    && rm -rf /var/lib/apt/lists/*

RUN conda config --add channels defaults \
  && conda config --add channels conda-forge \
  && conda config --add channels bioconda

RUN mamba install -c bioconda -y \
    blast=2.12.0 \
    && conda clean -afy

RUN mkdir /databases /input /output /work

# Import and format blast databases
WORKDIR /databases
RUN git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git
WORKDIR /databases/resfinder_db
RUN mkdir blast_db \
    && cat /databases/resfinder_db/*.fsa > /databases/resfinder_db/blast_db/concat_resfinder.fsa
WORKDIR /databases/resfinder_db/blast_db/
RUN makeblastdb -in concat_resfinder.fsa -dbtype nucl
ENV RESFINDER_DB="/databases/resfinder_db/blast_db/concat_resfinder.fsa"

WORKDIR /databases
RUN mkdir bacmet2_db \
    && cd bacmet2_db \
    && wget http://bacmet.biomedicine.gu.se/download/BacMet2_EXP_database.fasta \
    && dos2unix BacMet2_EXP_database.fasta \
    && sed -i '/^$/d' BacMet2_EXP_database.fasta
WORKDIR /databases/bacmet2_db
RUN makeblastdb -in BacMet2_EXP_database.fasta -dbtype prot
ENV BACMET2_DB="/databases/bacmet2_db/BacMet2_EXP_database.fasta"

WORKDIR /work

# for singularity compatibility
ENV LC_ALL=C
