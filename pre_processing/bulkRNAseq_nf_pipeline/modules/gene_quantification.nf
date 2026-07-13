process GENE_QUANT {
    tag "GENE_QUANT on $bam_files"
    publishDir "${params.out_dir}/03_gene_counts", pattern: "*.csv", mode: 'copy'
    label 'FEATURE_COUNTS'

    input:
    path bam_files

    output:
    path "${params.project_prefix}_gene_counts.csv", emit: quant_results
    
    script:
    """
    featureCounts -p -F GTF -s 0 -a $params.gene_gtf -o "${params.project_prefix}_gene_counts.csv" ${bam_files}
    """
}