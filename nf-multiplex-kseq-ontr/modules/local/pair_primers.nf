process PAIR_PRIMERS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'biocontainers/biopython:1.78' }"

    input:
    tuple val(meta), path(sites), path(genome)
    val min_amplicon
    val max_amplicon
    val min_amplicons_per_pair
    val max_amplicons_per_pair

    output:
    tuple val(meta), path("*.primer_pairs.txt"), emit: pairs
    path "*.pair_stats.txt"                     , emit: stats
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pair_primers.py \\
        --sites ${sites} \\
        --genome ${genome} \\
        --output ${prefix}.primer_pairs.txt \\
        --min-amplicon ${min_amplicon} \\
        --max-amplicon ${max_amplicon} \\
        --min-amplicons ${min_amplicons_per_pair} \\
        --max-amplicons ${max_amplicons_per_pair} \\
        --stats ${prefix}.pair_stats.txt \\
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
    touch ${prefix}.primer_pairs.txt
    touch ${prefix}.pair_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
