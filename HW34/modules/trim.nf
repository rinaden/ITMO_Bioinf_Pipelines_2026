process TRIM {

    tag "$id"

    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("*.fq.gz")

    script:
    if (reads.size() == 2) {
        """
        trim_galore --paired ${reads[0]} ${reads[1]}
        """
    } else {
        """
        trim_galore ${reads[0]}
        """
    }
}