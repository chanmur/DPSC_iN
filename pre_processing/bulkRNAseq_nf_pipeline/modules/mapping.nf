process UNIQUEMAPPING {
    tag "UniqueMapping $sample_id"
    label 'MAPPING'

    publishDir "${params.out_dir}/02_mapping", pattern: "${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_files)

    output:
    val sample_id, emit: sample_ids
    path "${sample_id}", emit: map_results
    tuple val(sample_id), val("${params.out_dir}/02_mapping/${sample_id}"), emit: map_results_path
    tuple val(sample_id), path("${sample_id}/*.bam"), emit: bam_files
    tuple val(sample_id), path("${sample_id}/*Log.final.out"), emit: log_final_files
    tuple val(sample_id), path("${sample_id}/*Log.out"), emit: log_out_files
    tuple val(sample_id), path("${sample_id}/*Log.progress.out"), emit: log_progress_files
    tuple val(sample_id), path("${sample_id}/*SJ.out.tab"), emit: sj_out_files
    tuple val(sample_id), path("${sample_id}/*_STARgenome"), emit: sample_genome_folders

    script:
    """
    mkdir ${sample_id}
    STAR --runThreadN 20 \
        --readFilesCommand gunzip -c \
        --outSAMattributes All \
        --outSAMtype BAM SortedByCoordinate \
        --genomeDir $params.star_index \
        --sjdbGTFfile $params.gene_gtf \
        --outFileNamePrefix "${sample_id}/${sample_id}_" \
        --outFilterMultimapNmax 1 \
        --outFilterMismatchNoverLmax 0.03  \
        --readFilesIn  $fastq_files
    """
}