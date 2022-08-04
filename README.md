[![Testing](https://github.com/EBI-Metagenomics/mobilomify/actions/workflows/test.yml/badge.svg)](https://github.com/EBI-Metagenomics/mobilomify/actions/workflows/test.yml)

# MoMofy
Mobilome Module for MGnify

Pipeline abstract.

## Setup

To install the pipeline we recommend Conda to manage the enviroment.

## Python dependencies

```python
$ pip install -r requirements-dev.txt
```

## Workflow

Pipeline diagram built with [cwlviewer](https://view.commonwl.org/) or similar.

# Tests

CWL tests are executed with [cwltest](https://github.com/common-workflow-language/cwltest).

Run:
```bash
$ cd tests
$ bash cwl-validation.sh
$ bash tools.sh
$ bash workflows.sh
```
