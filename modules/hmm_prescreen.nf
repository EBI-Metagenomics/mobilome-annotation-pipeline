process HMM_PRESCREEN {
    conda "bioconda::hmmer=3.3.2"
    
    input:
    tuple val(meta), path(faa_file)
    path hmm_profiles

    output:
    tuple val(meta), path("*_hmm_hits.txt"), emit: hmm_hits
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    hmmsearch \\
        --tblout ${meta.id}_hmm_hits.txt \\
        -E ${params.hmm_evalue ?: '1e-5'} \\
        --cpu ${task.cpus} \\
        ${args} \\
        ${hmm_profiles} \\
        ${faa_file}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmsearch -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//')
    END_VERSIONS
    """

}
