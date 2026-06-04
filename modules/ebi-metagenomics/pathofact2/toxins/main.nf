process PATHOFACT2_TOXINS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "microbiome-informatics/pathofact2_env:v1.0.4"

    input:
    tuple val(meta), path(fasta)
    path pathofact2_db

    output:
    tuple val(meta), path("*_classifier_toxins.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pathofact2'), eval("echo ${VERSION}"), topic: versions, emit: versions_pathofact2

    when:
    task.ext.when == null || task.ext.when

    script:
    VERSION = '1.0.4' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Handle compressed input
    def is_compressed = fasta.getExtension() == "gz"
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    def uncompress_input = is_compressed ? "gzip -c -d ${fasta} > ${fasta_name}" : ''

    // Handle compressed database
    def db_is_tarball = pathofact2_db.toString().endsWith('.tar.gz') || pathofact2_db.toString().endsWith('.tgz')
    def db_dir = db_is_tarball ? '.' : pathofact2_db
    def uncompress_db = db_is_tarball ? "tar -xzf ${pathofact2_db}" : ''

    """
    $uncompress_input
    $uncompress_db

    predict.py \\
        ${args} \\
        -s ${fasta_name} \\
        -m ${db_dir}/Models/TOX/final_model.joblib \\
        -v ${db_dir}/Models/TOX/ \\
        --cpus ${task.cpus} \\
        -o ${prefix}_classifier_toxins.tsv
    """

    stub:
    VERSION = '1.0.4' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    touch ${prefix}_classifier_toxins.tsv
    """
}
