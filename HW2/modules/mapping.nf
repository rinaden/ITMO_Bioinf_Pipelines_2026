process MAP {

    publishDir "${params.outdir}/mapping", mode: 'copy'

    input:
    tuple val(id), path(reads)
    path ref

    output:
    path "aligned.bam"

    script:
    if (reads.size() == 2) {
        """
        bwa index $ref
        bwa mem $ref ${reads[0]} ${reads[1]} | samtools sort -o aligned.bam
        samtools index aligned.bam
        """
    } else {
        """
        bwa index $ref
        bwa mem $ref ${reads[0]} | samtools sort -o aligned.bam
        samtools index aligned.bam
        """
    }
}
