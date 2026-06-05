process TRIM {

    tag "${meta.id}"

    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fq.gz")

    script:
    if(reads.size()==2){
        """
        trim_galore --paired ${reads[0]} ${reads[1]}
        """
    }
    else{
        """
        trim_galore ${reads[0]}
        """
    }
}