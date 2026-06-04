process INDEX_REFERENCE {

    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path fasta

    output:
    tuple path(fasta), path("${fasta}.fai")

    script:
    """
    samtools faidx $fasta
    """
}