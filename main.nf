#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
    nf-multiplex-kseq-ontr
========================================================================================
*/

include { GENERATE_PRIMERS        } from './modules/local/generate_primers'
include { FILTER_PRIMERS          } from './modules/local/filter_primers'
include { COUNT_SITES             } from './modules/local/count_sites'
include { PAIR_PRIMERS            } from './modules/local/pair_primers'
include { FIND_COMPATIBLE_SETS    } from './modules/local/find_compatible_sets'
include { OPTIMIZE_COVERAGE       } from './modules/local/optimize_coverage'
include { VALIDATE_PRIMERS        } from './modules/local/validate_primers'

workflow NFMULTIPLEXKSEQONTR {
    
    take:
    ch_genome
    
    main:
    GENERATE_PRIMERS(
        ch_genome,
        params.primer_length
    )
    
    FILTER_PRIMERS(
        GENERATE_PRIMERS.out.primers,
        params.min_tm,
        params.max_tm,
        params.min_gc,
        params.max_gc
    )
    
    COUNT_SITES(
        FILTER_PRIMERS.out.filtered.join(ch_genome),
        params.min_binding_sites,
        params.max_binding_sites
    )
    
    PAIR_PRIMERS(
        COUNT_SITES.out.sites.join(ch_genome),
        params.min_amplicon_size,
        params.max_amplicon_size,
        params.min_amplicons_per_pair,
        params.max_amplicons_per_pair
    )
    
    FIND_COMPATIBLE_SETS(
        PAIR_PRIMERS.out.pairs,
        params.min_primer_pairs,
        params.max_primer_pairs
    )
    
    OPTIMIZE_COVERAGE(
        FIND_COMPATIBLE_SETS.out.sets,
        params.target_snps,
        params.min_snps
    )
    
    VALIDATE_PRIMERS(
        OPTIMIZE_COVERAGE.out.optimized.join(ch_genome)
    )
    
    emit:
    final_primers = VALIDATE_PRIMERS.out.final_primers
    validation_report = VALIDATE_PRIMERS.out.report
}

workflow {
    
    main:
    if (!params.genome) {
        error "Please provide a genome file with --genome"
    }
    
    def meta = [id: params.sample_id ?: 'sample']
    def ch_genome = channel.fromPath(params.genome, checkIfExists: true)
                           .map { file -> [meta, file] }
    
    NFMULTIPLEXKSEQONTR(ch_genome)
    
    workflow.onComplete {
        def msg = """\
            Pipeline execution summary
            ---------------------------
            Completed at: ${workflow.complete}
            Duration    : ${workflow.duration}
            Success     : ${workflow.success}
            workDir     : ${workflow.workDir}
            exit status : ${workflow.exitStatus}
            """.stripIndent()
        
        println msg
        
        if (workflow.success) {
            println "Results are in: ${params.outdir}"
        }
    }
    
    workflow.onError {
        println "Error: ${workflow.errorMessage}"
    }
}
