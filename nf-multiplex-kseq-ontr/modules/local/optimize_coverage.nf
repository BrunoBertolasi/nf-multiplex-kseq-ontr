process OPTIMIZE_COVERAGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'biocontainers/biopython:1.78' }"

    input:
    tuple val(meta), path(sets)
    val target_snps
    val min_snps

    output:
    tuple val(meta), path("*.optimized_set.txt"), emit: optimized
    path "*.optimization_report.txt"             , emit: report
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    optimize_coverage.py \\
        --sets ${sets} \\
        --output ${prefix}.optimized_set.txt \\
        --target-snps ${target_snps} \\
        --min-snps ${min_snps} \\
        --report ${prefix}.optimization_report.txt \\
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
    touch ${prefix}.optimized_set.txt
    touch ${prefix}.optimization_report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
