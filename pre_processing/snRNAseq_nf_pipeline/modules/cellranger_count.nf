process COUNT {
    tag "COUNT on ${sample_id}"
    
    publishDir "${params.out_dir}/01_count_output/", pattern: "${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_path)

    output:
    val sample_id, emit: sample_ids
    path "${sample_id}", emit: count_results
    tuple val(sample_id), path("${sample_id}"), emit: count_dir_tuples
    tuple val(sample_id), path("${sample_id}/outs/possorted_genome_bam.bam"), emit: bam_files
    tuple val(sample_id), path("${sample_id}/outs/possorted_genome_bam.bam.bai"), emit: bai_files
    tuple val(sample_id), path("${sample_id}/outs/filtered_feature_bc_matrix"), emit: filtered_matrices
    tuple val(sample_id), path("${sample_id}/outs/raw_feature_bc_matrix"), emit: raw_matrices
    tuple val(sample_id), path("${sample_id}/outs/web_summary.html"), emit: web_summary
    tuple val(sample_id), path("${sample_id}/outs/metrics_summary.csv"), emit: metrics_summary

    script:
    """
    mkdir fastq_dir
    rsync -azvhPt ${fastq_path.join(' ')} fastq_dir/ 
    
    cellranger count --id=${sample_id} \
           --transcriptome=${params.transcriptome_index}\
           --fastqs=fastq_dir \
           --sample=${sample_id} \
           --create-bam=true
    """
}


process ORGANIZE_COUNT_OUTPUT {
    tag "ORGANIZE_COUNT_OUTPUT on ${sample_id}"
    
    publishDir "${params.out_dir}/01_count_output/organized/", mode: 'copy'

    input:
    tuple val(sample_id), path(count_dir)

    output:
    val sample_id, emit: sample_ids
    path "filtered_matrices/${sample_id}", emit: filtered_matrices
    path "bam_files/${sample_id}_possorted_genome_bam.bam", emit:bam_files
    path "bam_files/${sample_id}_possorted_genome_bam.bam.bai", emit: bai_files
    path "metrics/${sample_id}_metrics_summary.csv", emit: metrics_summary
    path "metrics/${sample_id}_web_summary.html", emit: web_summary

    script:
    """
    mkdir -p filtered_matrices/${sample_id} bam_files metrics

    # Matric files
    if [ -d "${count_dir}/outs/filtered_feature_bc_matrix" ]; then
        rsync -azvhP ${count_dir}/outs/filtered_feature_bc_matrix/ filtered_matrices/${sample_id}/
    fi

    # Bam and index files
    if [ -f "${count_dir}/outs/possorted_genome_bam.bam" ]; then
        rsync -azvhP ${count_dir}/outs/possorted_genome_bam.bam bam_files/${sample_id}_possorted_genome_bam.bam
    fi
    if [ -f "${count_dir}/outs/possorted_genome_bam.bam.bai" ]; then
        rsync -azvhP ${count_dir}/outs/possorted_genome_bam.bam.bai bam_files/${sample_id}_possorted_genome_bam.bam.bai
    fi

    # Metrics and web summary files
    if [ -f "${count_dir}/outs/metrics_summary.csv" ]; then
        rsync -azvhP ${count_dir}/outs/metrics_summary.csv metrics/${sample_id}_metrics_summary.csv
    fi
    if [ -f "${count_dir}/outs/web_summary.html" ]; then
        rsync -azvhP ${count_dir}/outs/web_summary.html metrics/${sample_id}_web_summary.html
    fi
    """
}