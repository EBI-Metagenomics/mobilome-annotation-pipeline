#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process subfastq {
    publishDir "$launchDir/preprocessing/subreads", mode: 'copy'

    memory "8 GB"
    cpus 1

    input:
      path fq_1, name: 'input_1.fastq.gz'
      path fq_2, name: 'input_2.fastq.gz'

    output:
      path 'sub_1.fq.gz', emit: pair_1
      path 'sub_2.fq.gz', emit: pair_2

    script:
      """
      #!/bin/bash
      zcat input_1.fastq.gz | head -40000000 > sub_1.fq && gzip sub_1.fq
      zcat input_2.fastq.gz | head -40000000 > sub_2.fq && gzip sub_2.fq
      """
}
