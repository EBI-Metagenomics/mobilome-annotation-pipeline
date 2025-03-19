process ISESCAN {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/microbiome-informatics/isescan-v1.7.2.3'

    input:
    tuple val(meta), path(assembly_file)

    output:
    tuple val(meta), path('isescan_results/1kb_contigs.fasta.tsv'), emit: iss_tsv

    script:
    """
    isescan.py --seqfile ${assembly_file} \\
    --output isescan_results \\
    --nthread ${task.cpus}

    if ls -l isescan_results/1kb_contigs.fasta.tsv 2>/dev/null | grep -q .
    then
        echo 'ISEScan results exists'
    else
    echo 'ISEScan found 0 insertion sequences in input file... generating dummy files'
        touch isescan_results/1kb_contigs.fasta.tsv
    fi
    """
}
