process INDEX_BAM {
    tag "INDEX_BAM on $bam_file"
    publishDir "$bam_path", pattern: "*.bai", mode: 'copy'

    input:
    tuple val(sample_id), val(bam_path), path(bam_file)

    output:
    val sample_id, emit: sample_ids
    tuple val(sample_id), val(bam_path), emit: map_results_path
    tuple val(sample_id), path(bam_file), emit: bam_files
    tuple val(sample_id), path("*.bai"), emit: bai_files

    script:
    """
    samtools index -b $bam_file
    """
}

process GENERATE_BIGWIG {
    tag "GENERATE_BIGWIG on $bam_file and $bai_file"
    label 'BIGWIG'
    publishDir "$bam_path", pattern: "*.bw", mode: 'copy'

    input:
    tuple val(sample_id), val(bam_path), path(bam_file), path(bai_file)
    
    output:
    val sample_id, emit: sample_ids
    path "*.bw", emit: bigwig_files
    tuple val(sample_id), val(bam_path), emit: map_results_path
    tuple val(sample_id), path("*.bw"), emit: bai_files

    script:
    """
    bamCoverage --normalizeUsing RPKM -b $bam_file -o "${sample_id}.bw"
    """
}