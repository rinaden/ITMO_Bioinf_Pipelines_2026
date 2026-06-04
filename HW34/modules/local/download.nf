process DOWNLOAD {

    tag "$id"

    publishDir "${params.outdir}/download", mode: 'copy'

    input:
    val id

    output:
    path "*.fastq.gz"

    script:
    """
    fasterq-dump $id --split-files -e 4
    gzip *.fastq
    """
}