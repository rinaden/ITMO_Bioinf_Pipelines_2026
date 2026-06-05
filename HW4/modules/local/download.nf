process DOWNLOAD {

    tag "${meta.id}"

    publishDir "${params.outdir}/download", mode: 'copy'

    input:
    tuple val(meta), val(sra)

    output:
    tuple val(meta), path("*.fastq.gz")

    script:
    """
    fasterq-dump ${sra} --split-files -e ${task.cpus}
    gzip *.fastq
    """
}