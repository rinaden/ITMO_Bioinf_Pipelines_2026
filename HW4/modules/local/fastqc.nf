process FASTQC {

    publishDir "${params.outdir}/fastqc", mode: 'copy'

    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_fastqc.zip")

    script:
    """
    fastqc ${reads.join(' ')}
    """
}