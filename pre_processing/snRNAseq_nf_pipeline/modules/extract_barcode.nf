process EXTRACT_BARCODE {
    tag "BARCODE on ${sample_id}"
    
    publishDir "${params.out_dir}/02_barcode_output/", pattern: "${sample_id}_*", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_path)

    output:
    val sample_id, emit: sample_ids
    tuple val(sample_id), path("${sample_id}_barcodes.tsv"), emit: barcodes
    tuple val(sample_id), path("${sample_id}_non_libreads.fa"), emit: nonlib

    script:
    """
    extract_barcodes ${sample_id} ${fastq_path.join(' ')}
    """
}