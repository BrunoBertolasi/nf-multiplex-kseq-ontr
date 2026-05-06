process FILTER_PRIMERS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'biocontainers/biopython:1.78' }"

    input:
    tuple val(meta), path(primers)
    val min_tm
    val max_tm
    val min_gc
    val max_gc

    output:
    tuple val(meta), path("*.filtered_primers.txt"), emit: filtered
    path "*.filter_stats.txt"                       , emit: stats
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    filter_properties.py \\
        --input ${primers} \\
        --output ${prefix}.filtered_primers.txt \\
        --min-tm ${min_tm} \\
        --max-tm ${max_tm} \\
        --min-gc ${min_gc} \\
        --max-gc ${max_gc} \\
        --stats ${prefix}.filter_stats.txt \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filtered_primers.txt
    touch ${prefix}.filter_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
