process COUNT_SITES {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'biocontainers/biopython:1.78' }"

    input:
    tuple val(meta), path(primers), path(genome)
    val min_sites
    val max_sites

    output:
    tuple val(meta), path("*.binding_sites.txt"), emit: sites
    path "*.site_stats.txt"                      , emit: stats
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    count_sites.py \\
        --primers ${primers} \\
        --genome ${genome} \\
        --output ${prefix}.binding_sites.txt \\
        --min-sites ${min_sites} \\
        --max-sites ${max_sites} \\
        --stats ${prefix}.site_stats.txt \\
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
    touch ${prefix}.binding_sites.txt
    touch ${prefix}.site_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
