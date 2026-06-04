process ASSEMBLE {

    publishDir "${params.outdir}/assembly", mode: 'copy'

    input:
    tuple val(id), path(reads)

    output:
    path "contigs.fasta"

    script:
    if (reads.size() == 2) {
        """
        spades.py -1 ${reads[0]} -2 ${reads[1]} -o spades_out
        cp spades_out/contigs.fasta .
        """
    } else {
        """
        spades.py -s ${reads[0]} -o spades_out
        cp spades_out/contigs.fasta .
        """
    }
}
