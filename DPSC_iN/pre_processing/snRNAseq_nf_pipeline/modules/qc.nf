process MULTIQC {
    tag "MULTIQC"

    publishDir "${params.out_dir}/qc/multi_qc", mode: 'copy', overwrite: true

    input:
    path qc_files

    output:
    path "multiqc_data"
    path "multiqc_report.html"

    script:
    """
    multiqc -d -s -f ${qc_files}
    """
}