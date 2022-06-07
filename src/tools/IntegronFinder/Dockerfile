# Dockerfile for integronfinder

FROM quay.io/microbiome-informatics/miniconda:4.11.0

ENV VERSION="2.0.2"

LABEL maintainer="Microbiome Informatics Team www.ebi.ac.uk/metagenomics"
LABEL software="integronfinder"
LABEL version="${VERSION}"
LABEL description="Tool to find integrons in DNA sequences (contigs)"
LABEL website="https://github.com/gem-pasteur/Integron_Finder"
LABEL license="GPLv3.0"


# Setting up channels
RUN conda config --add channels defaults \
  && conda config --add channels conda-forge \
  && conda config --add channels bioconda

# install dependencies
RUN conda install -c bioconda -y infernal=1.1.4 hmmer=3.3.2 prodigal=2.6.3 \
   && conda clean -afy

# install integronfinder via pip
RUN pip install integron_finder==${VERSION}

# for singularity compatibility
ENV LC_ALL=C

# setting workdir
RUN mkdir /input /output /work
WORKDIR /work

# simple executor
CMD ["/bin/bash", "-c"]

