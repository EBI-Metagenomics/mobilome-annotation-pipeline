process INTERPROSCAN {
    tag "$meta.id"
    label 'process_long'

    container 'microbiome-informatics/interproscan:5.76-107.0_patch1'

    containerOptions {
        def containerArgs = []
        def mountArg = (workflow.containerEngine == 'singularity') ? "--bind" : "--volume"

        containerArgs << "${mountArg} ${task.workDir}/${interproscan_db}/data:/opt/interproscan/data"

        if ( params.interpro_licensed_software ) {
            def licensedSoftwarePath = "${task.workDir}/${interproscan_db}/licensed"
            containerArgs << "${mountArg} ${licensedSoftwarePath}:/opt/interproscan/licensed"
            // This override is needed otherwise it fails because this path seems to be hardcoded in the container
            containerArgs << "${mountArg} ${licensedSoftwarePath}/signalp:/usr/opt/www/pub/CBS/services/SignalP-4.1/signalp-4.1"
        }

        return containerArgs.join(' ')
    }

    input:
    tuple val(meta), path(fasta)
    tuple path(interproscan_db), val(db_version)

    output:
    tuple val(meta), path('*.tsv.gz') , emit: tsv
    tuple val(meta), path('*.xml.gz') , emit: xml
    tuple val(meta), path('*.gff3.gz'), emit: gff3
    tuple val(meta), path('*.json.gz'), emit: json
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.extension == "gz"
    def fasta_file_name = fasta.name - ~/\.gz$/
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_file_name}
    fi

    # Set the max memory for the JVM
    export JAVA_OPTS="-Xmx${task.memory.toGiga()}G"

    # -dp (disable precalculation) is on so no online dependency
    interproscan.sh \\
        -cpu $task.cpus \\
        -i ${fasta_file_name} \\
        -dp \\
        ${args} \\
        --output-file-base ${prefix}

    gzip ${prefix}.{tsv,xml,gff3,json}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        InterProScan: \$(interproscan.sh --version | grep -o "InterProScan version [0-9.-]*" | sed "s/InterProScan version //")
        InterProScan database: $db_version
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' > ${prefix}.{tsv,xml,gff3,json}
    gzip ${prefix}.{tsv,xml,gff3,json}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        InterProScan: \$(interproscan.sh --version | grep -o "InterProScan version [0-9.-]*" | sed "s/InterProScan version //")
        InterProScan database: $db_version
    END_VERSIONS
    """
}
