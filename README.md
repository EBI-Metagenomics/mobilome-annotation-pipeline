[![Testing](https://github.com/EBI-Metagenomics/mobilomify/actions/workflows/test.yml/badge.svg)](https://github.com/EBI-Metagenomics/mobilomify/actions/workflows/test.yml)

# MoMofy
Mobilome Module for MGnify

MoMofy is a wraper that ibtegrates the ouptput of differeten tools designed fot the prediction of integrative mobile genetic elements in prokaryotic genomes and metagenomes. If you use MoMofy, please cite the papers of the original tools:

1) ISEScan v1.7.2.3 [Xie et al., Bioinformatics, 2017](https://doi.org/10.1093/bioinformatics/btx433)
2) PaliDIS v3.1.2 [Carr et al., biorxiv, 2022](https://doi.org/10.1101/2022.06.27.497710)
3) IntegronFinder2 v2.0.2 [Néron et al., Microorganisms, 2022](https://doi.org/10.3390/microorganisms10040700)
4) ICEfinder v1.0 [Liu et al., Nucleic Acids Research, 2019](https://doi.org/10.1093/nar/gky1123)

Databases:
- MobileOG-DB Beatrix 1.6 v1 [Brown et al., Appl Environ Microbiol, 2022](https://doi.org/10.1128/aem.00991-22)


## Setup

This workflow was implemented to be ran through Nextflow v21.10.0. Most of the tools are available on quay.io and no install is needed. 

In the case of ICEfinder, the user will need to contact the author to get their own copy of the software, visit the [ICEfinder website](https://bioinfo-mml.sjtu.edu.cn/ICEfinder/download.html). Once this is done, use the Dockerfile template provided in this repo (templates/icefinder/Dockerfile) to built your own container and set up the corresponding parameter on the nextflow.config file (line XX).


## Workflow

|workflow|\


# Tests

Nextflow tests are executed with [nf-test](https://github.com/askimed/nf-test).

Run:
```bash
$ cd test
$ nf-test test *.nf.test
```


.. |workflow| image:: media/momofy_workflow.png
