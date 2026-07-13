process FASTQC {
    tag "FASTQC on $sample_id"
    
    publishDir "${params.out_dir}/01_qc/fastqc", pattern: "${sample_id}_fastqc_logs", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_files)

    output:
    val sample_id, emit: sample_ids
    path "${sample_id}_fastqc_logs", emit: fastqc_out_files
    path "${sample_id}_fastqc_logs/*.html", emit: html_files
    path "${sample_id}_fastqc_logs/*.zip", emit: zip_files 

    script:
    """
    mkdir ${sample_id}_fastqc_logs
    fastqc -o "${sample_id}_fastqc_logs" -f fastq -q ${fastq_files}
    """
}

process RSEQC_INFER_STRAND {
    tag "RSEQC_INFER_STRAND on $sample_id"
    
    publishDir "${params.out_dir}/01_qc/reseqc", pattern: "${sample_id}_infer_experiment.txt", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    val sample_id, emit: sample_ids
    path "${sample_id}_infer_experiment.txt", emit: reseq_file 

    script:
    """
    infer_experiment.py -i $bam -r $params.rseqc_bed > ${sample_id}_infer_experiment.txt
    """
}

process MULTIQC {
    tag "MULTIQC"

    publishDir "${params.out_dir}/01_qc/multi_qc", mode: 'copy', overwrite: true

    input:
    path "*"

    output:
    path "multiqc_data"
    path "multiqc_report.html"

    script:
    """
    multiqc -d -s -f .
    """
}