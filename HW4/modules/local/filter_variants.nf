process FILTER_VARIANTS {

    tag "${meta.id}"

    publishDir "${params.outdir}/filtered", mode: 'copy'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.filtered.vcf.gz")

    script:
    """
    bcftools filter \
        -i 'QUAL>20' \
        $vcf \
        -Oz \
        -o ${meta.id}.filtered.vcf.gz

    tabix -p vcf ${meta.id}.filtered.vcf.gz
    """

    stub:
    """
    touch ${meta.id}.filtered.vcf.gz
    touch ${meta.id}.filtered.vcf.gz.tbi
    """
}