process CALLING {

    tag "${sample_id}"

    publishDir "${params.outdir}/calling", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz")

    script:
    """
    samtools index $bam

    bcftools mpileup \
        -f $reference \
        $bam \
    | bcftools call \
        -mv -Oz \
        -o ${sample_id}.vcf.gz

    tabix -p vcf ${sample_id}.vcf.gz
    """
}