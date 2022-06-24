# Dockerfile for blast
FROM quay.io/microbiome-informatics/miniconda:4.11.0


ENV VERSION="2.12.0"

LABEL maintainer="Microbiome Informatics Team www.ebi.ac.uk/metagenomics"
LABEL software="Blast"
LABEL version="${VERSION}"
LABEL description="Blast tool"
LABEL webpage="https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download"
LABEL license="https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/scripts/projects/blast/LICENSE"

# Dependencies
RUN conda config --add channels bioconda

RUN mamba install -c bioconda -y \
    blast=2.12.0 \
    && conda clean -afy

# for singularity compatibility
ENV LC_ALL=C
