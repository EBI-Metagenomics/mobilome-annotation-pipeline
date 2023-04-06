#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process database_format {
    publishDir "$projectDir/databases", mode: 'copy'
    stageInMode = 'copy'

    memory "8 GB"
    cpus 4

    container "quay.io/biocontainers/diamond:2.0.12--hdcc8f71_0"

    input:
        path mobileog

    output:
        path 'mobileOG_beatrix1.6.dmnd'

    script:
        """
	diamond makedb --in ${mobileog} --db mobileOG_beatrix1.6.dmnd
        """
}
workflow {
    Channel.fromPath('mobileOG-db_beatrix-1.6.All.faa') | database_format
}
