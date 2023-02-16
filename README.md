[![Testing](https://github.com/EBI-Metagenomics/mobilomify/actions/workflows/test.yml/badge.svg)](https://github.com/EBI-Metagenomics/mobilomify/actions/workflows/test.yml)

# MoMofy
Mobilome Module for MGnify

MoMofy is a wraper that ibtegrates the ouptput of differeten tools designed fot the prediction of integrative mobile genetic elements in prokaryotic genomes and metagenomes. If you use MoMofy, please cite the papers of the original tools:

1) ISEScan v1.7.2.3
2) PaliDIS v3.1.2
3) IntegronFinder2 v2.0.2
4) ICEfinder v1.0

Databases:
- MobileOG-DB


## Setup

This workflow was implemented to run through Nextflow v21.10.0. Most of the tools are available on quay.io and no install is needed. In the case of ICEfinder, the user will need to contact the author to get their own copy of the software, visit the [ICEfinder website](https://bioinfo-mml.sjtu.edu.cn/ICEfinder/download.html). Once it is done, use the Dockerfile template provided in this repo (templates/icefinder/Dockerfile) to built your own container and set up the corresponding parameter on the nextflow.config file (line XX).



## Python dependencies

```python
$ pip install -r requirements-dev.txt
```

## Workflow



# Tests

Nextflow tests are executed with [nf-test](https://github.com/askimed/nf-test).

Run:
```bash
$ cd test
$ nf-test test *.nf.test
```
