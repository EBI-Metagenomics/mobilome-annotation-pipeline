

# MoMofy
Mobilome Module for MGnify

Introductory sentence...

MoMofy is a wraper that integrates the ouptput of differeten tools designed for the prediction of integrative mobile genetic elements in prokaryotic genomes and metagenomes. 


## Contents
- [ Workflow ](#wf)
- [ Setup ](#sp)
- [ MoMofy install and dependencies ](#install)
- [ Inputs ](#in)
- [ Outputs ](#out)
- [ Tests ](#test)
- [ Citation ](#cite)


<a name="wf"></a>
## Workflow

<img src="media/momofy_workflow.png" width="1600"/>


<a name="sp"></a>
## Setup

This workflow is implemented to be run through [Nextflow](https://www.nextflow.io/) and use [Docker](https://www.docker.com/) to build the image of ICEfinder. 

- Install [Nextflow version >=21.10](https://www.nextflow.io/docs/latest/getstarted.html#installation)
- Install [Docker](https://docs.docker.com/get-docker/)


<a name="install"></a>
## MoMofy install and dependencies

To instal MoMofy, clone this repo by:

```bash
$ git clone https://github.com/EBI-Metagenomics/momofy.git
$ cd momofy
$ bash setup.sh
```

Most of the tools are available on [quay.io](https://quay.io) and no install is needed. 

In the case of ICEfinder, the user will need to contact the author to get their own copy of the software, visit the [ICEfinder website](https://bioinfo-mml.sjtu.edu.cn/ICEfinder/download.html). Once you have the ICEfinder_linux.tar.gz tarball, move it to momofy/templates/icefinder/ and build the docker image:

```bash
$ mv ICEfinder_linux.tar.gz /PATH/momofy/templates/icefinder/
$ cd /PATH/momofy/templates/icefinder/
$ docker build -t my_icefinder .
```

<a name="in"></a>
## Inputs

To run MoMofy create a directory per sample and launch the tool from the sample directory. If you have many samples, you can use a list of sample IDs and iterate on it to create multiple directories and subdirectories, prepare your data inputs and then run MoMofy. Here an example on how to create softlinks to the assembly files:

```bash
$ for sample in $(cat samples.list); do (mkdir -p $sample/$raw_data && cd $sample/$raw_data && ln -s /path/to/assemblies/$sample.fasta $sample/$raw_data/contigs.fasta . ); done
```

Two input files are mandatory in a folder called `raw_data`:
- (meta)genomic assembly file in fasta format (uncompress)
- proteins predicted in fasta format (uncompress)

Optional input files to run PaliDIS on metagenomic assemblies (also in `raw_data` folder):
- read 1 fastq file (gzip compress)
- read 2 fastq file (gzip compress)

The inputs can be symbolic links and your `raw_data` folder should look like this:

```bash
$ tree raw_data
raw_data/
├── contigs.fasta
├── file_1.fastq.gz
├── file_2.fastq.gz
└── proteins.faa
```

Basic usage:

```bash
$ cd 
$ nextflow run /PATH/momofy/momofy.nf --assembly raw_data/contigs.fasta --cds_annot raw_data/proteins.faa -with-docker my_icefinder
```

Running PaliDIS requires as input the Illumina paired-end reads used to generate the metagenomic assembly. By default, this option is turned off. To run PaliDIS add the following parameters to the command:

```bash
$ cd /PATH/momofy
$ nextflow run momofy.nf --assembly raw_data/contigs.fasta --cds_annot raw_data/proteins.faa --palidis true --read_1 raw_data/file_1.fastq.gz --read_2 raw_data/file_2.fastq.gz -with-docker my_icefinder
```

<a name="out"></a>
## Outputs

The main output directory of MoMofy is `MoMofy_results` and contain four files:

```bash
$ tree MoMofy_results
MoMofy_results/
├── discarded_iss.txt
├── momofy_predictions.fna
├── momofy_predictions.gff
└── nested_integrons.txt
```

Additionally, you will see the directories containing the main outputs of each tool. This is a minimal example omiting the input directory:

```bash
$ tree ./
├── icefinder_results
│   ├── gbk
│   │   └── contig.prokka_1.gbk
│   ├── input.list
│   ├── result
│   │   ├── icf_concat.fasta
│   │   └── icf_concat.summary
│   └── tmp
│       ├── contig.prokka_1
├── integron_results
│   └── Results_Integron_Finder_contigs
│       ├── contig_1.gbk
│       └── contigs.summary
├── isescan_results
│   ├── contigs.fasta.is.fna
│   └── contigs.fasta.tsv
├── mobileog_results
│   └── blastp_out.tsv
└── preprocessing
    ├── 1kb_contigs.fasta
    ├── 5kb_contigs.fasta
    ├── contigID.map
    ├── prokka_out
    │   └── contigs.gbk
    └── subreads
        ├── sub_1.fq.gz
        └── sub_2.fq.gz

```

<a name="test"></a>
## Tests

Nextflow tests are executed with [nf-test](https://github.com/askimed/nf-test).

Run:
```bash
$ cd test
$ nf-test test *.nf.test
```

<a name="cite"></a>
## Citation

If you use MoMofy on your data analysis, please cite:

XXXXX


MoMofy is a wrapper that integrates the output of the following tools and DBs:

1) ISEScan v1.7.2.3 [Xie et al., Bioinformatics, 2017](https://doi.org/10.1093/bioinformatics/btx433)
2) PaliDIS v3.1.2 [Carr et al., biorxiv, 2022](https://doi.org/10.1101/2022.06.27.497710)
3) IntegronFinder2 v2.0.2 [Néron et al., Microorganisms, 2022](https://doi.org/10.3390/microorganisms10040700)
4) ICEfinder v1.0 [Liu et al., Nucleic Acids Research, 2019](https://doi.org/10.1093/nar/gky1123)

Databases:
- MobileOG-DB Beatrix 1.6 v1 [Brown et al., Appl Environ Microbiol, 2022](https://doi.org/10.1128/aem.00991-22)
