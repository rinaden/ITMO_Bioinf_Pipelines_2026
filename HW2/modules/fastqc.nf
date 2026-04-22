process FASTQC {

    publishDir "${params.outdir}/fastqc", mode: 'copy'

    tag "$reads"

    input:
    tuple val(id), path(reads)

    output:
    path "*_fastqc.zip"

    script:
    """
    fastqc ${reads.join(' ')}
    """
}
