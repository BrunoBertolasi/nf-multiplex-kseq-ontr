process FIND_COMPATIBLE_SETS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'biocontainers/biopython:1.78' }"

    input:
    tuple val(meta), path(pairs)
    val min_pairs
    val max_pairs

    output:
    tuple val(meta), path("*.compatible_sets.txt"), emit: sets
    path "*.compatibility_stats.txt"               , emit: stats
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    find_compatible_sets.py \\
        --pairs ${pairs} \\
        --output ${prefix}.compatible_sets.txt \\
        --min-pairs ${min_pairs} \\
        --max-pairs ${max_pairs} \\
        --stats ${prefix}.compatibility_stats.txt \\
        --threads ${task.cpus} \\
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
    touch ${prefix}.compatible_sets.txt
    touch ${prefix}.compatibility_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
