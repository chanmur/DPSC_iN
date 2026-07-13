#!/usr/bin/env nextflow

nextflow.enable.dsl=2


// importing modules
include {FASTQC; MULTIQC; RSEQC_INFER_STRAND} from './modules/quality_control'
include {UNIQUEMAPPING} from './modules/mapping'
include {INDEX_BAM; GENERATE_BIGWIG} from './modules/bam_processing'
include {GENE_QUANT} from './modules/gene_quantification'

workflow qc_mapping {
    take:
    fastq_read_pairs

    main:
    fastq_read_pairs | (FASTQC & UNIQUEMAPPING)

    multi_qc_ch = FASTQC.out.fastqc_out_files
                    .mix(UNIQUEMAPPING.out.map_results)

    multi_qc_ch | collect | MULTIQC

    emit:
    map_path_ch = UNIQUEMAPPING.out.map_results_path
    map_bam_ch = UNIQUEMAPPING.out.bam_files
}

workflow infer_strandedness {
    take:
    map_bam_ch

    main:
    map_bam_ch | RSEQC_INFER_STRAND

}

workflow index_bigwig {
    take:
    map_path_ch
    map_bam_ch

    main:
    map_ch = map_path_ch
                .join(map_bam_ch, remainder: false)

    map_ch | INDEX_BAM

    index_ch = INDEX_BAM.out.map_results_path
                .join(INDEX_BAM.out.bam_files, remainder: false)
                .join(INDEX_BAM.out.bai_files, remainder: false)
    
    index_ch | GENERATE_BIGWIG
}

workflow quant {
    take:
    map_bam_ch

    main:
    bam_files = map_bam_ch
        .map { it -> it[1] }
        .collect() 
        
    bam_files | GENE_QUANT

}

workflow preprocessing {
    take:
    fastq_read_pairs

    main:
    qc_mapping_out = qc_mapping(fastq_read_pairs)

    map_path_ch = qc_mapping_out.map_path_ch
    map_bam_ch = qc_mapping_out.map_bam_ch

    infer_strandedness(map_bam_ch)
    index_bigwig(map_path_ch, map_bam_ch)
    quant(map_bam_ch)
}

// Main workflow
workflow {
    fastq_read_pairs = channel.fromFilePairs(
        params.fastq, 
        checkIfExists: true)
        .map {sample_id, read_pair -> 
            def parts = sample_id.tokenize('_')
            def donor = parts[2]
            def cell_line = parts[4]
            def new_id = "${donor}_${cell_line}"
            tuple(new_id, read_pair)
        }

    preprocessing(fastq_read_pairs)
}