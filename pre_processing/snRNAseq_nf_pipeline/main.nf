#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// importing modules
include {COUNT; ORGANIZE_COUNT_OUTPUT} from './modules/cellranger_count'
include {MULTIQC} from './modules/qc'
include {EXTRACT_BARCODE} from './modules/extract_barcode'

// The count workflow below assume that all the fastq files from across the samples are all in the same folder. If they are not, it needs to be adjusted accordingly

workflow count_snrna {
    take:
    sn_rnaseq_fastq

    main:
    
    sn_rnaseq_fastq | COUNT

    COUNT.out.count_dir_tuples | ORGANIZE_COUNT_OUTPUT

    ORGANIZE_COUNT_OUTPUT.out.metrics_summary 
    .collect()
    | MULTIQC
}

workflow lineage_barcodes {
    take: 
    bc_fastq

    main:
    bc_fastq | EXTRACT_BARCODE
}

workflow {
    channel.fromFilePairs(
        params.fastq, 
        checkIfExists: true, size: -1) { file -> 
            def name = file.getName()
            def match = name =~ /(.*)_S\d+_L\d+_R\d+_001\.fastq\.gz/
            def sample_id = match ? match[0][1] : { throw new IllegalArgumentException("Filename does not match expected pattern: $name") }()
        }
        .branch { sample_id, files ->
            bc_fastq: sample_id.endsWith('_barcode')
                return tuple(sample_id.replaceAll(/_barcode$/, ''), files)
            
            sn_rnaseq_fastq: !sample_id.endsWith('_barcode')
                return tuple(sample_id, files)
        }
        .set { fastq_read_pairs }


    count_snrna(fastq_read_pairs.sn_rnaseq_fastq)
    //lineage_barcodes(fastq_read_pairs.bc_fastq)
}
