process SANNTIS {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sanntis:0.9.4.1--pyhdfd78af_0'
        : 'biocontainers/sanntis:0.9.4.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(interproscan), path(gbk), path(faa)

    output:
    tuple val(meta), path("*_sanntis.gff.gz"), emit: gff
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (faa && gbk) {
        error("SanntiS supports either a GBK or a FAA as the secondary input.")
    }
    if (!gbk && !faa) {
        error("SanntiS needs either a GBK or FAA as input.")
    }

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_ips_compressed = interproscan.name.endsWith(".gz")
    def  interproscan_file = interproscan ? interproscan.name.replace(".gz", "") : is_ips_compressed

    // Handle the GBK or FAA
    def input_arg = ""
    def input_file = ""
    def decompress = ""
    def clean_decompressed = ""
    if (gbk) {
        input_file = gbk.name
        if (gbk.name.endsWith(".gz")) {
            input_file = input_file.replace(".gz", "")
            decompress = "gzip -c -d ${gbk} | grep -v \"/protein_id=\" > ${input_file}"
            clean_decompressed = "rm -f ${input_file}"
        }
        input_arg = "${input_file}"
    }
    if (faa) {
        input_file = faa.name
        if (faa.name.endsWith(".gz")) {
            input_file = input_file.replace(".gz", "")
            decompress = "gzip -c -d ${faa} > ${input_file}"
            clean_decompressed = "rm -f ${input_file}"
        }
        input_arg = "--is_protein ${input_file}"
    }
    """
    if [ "${is_ips_compressed}" == "true" ]; then
         gzip -c -d ${interproscan} > ${interproscan_file}
    fi

    ${decompress}

    sanntis \\
        --ip-file ${interproscan_file} \\
        --outfile ${prefix}_sanntis.gff \\
        --cpu ${task.cpus} \\
        ${args} \\
        ${input_arg}

    gzip ${prefix}_sanntis.gff

    # Purge the decompressed files to save storage
    ${clean_decompressed}
    if [ "${is_ips_compressed}" == "true" ]; then
         rm -f ${interproscan_file}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SanntiS: \$(sanntis --version | sed "s/SanntiS\\ //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_sanntis.gff
    gzip ${prefix}_sanntis.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SanntiS: \$(sanntis --version | sed "s/SanntiS\\ //g")
    END_VERSIONS
    """
}
