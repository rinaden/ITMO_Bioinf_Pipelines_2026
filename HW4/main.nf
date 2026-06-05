nextflow.enable.dsl=2

include { TRIM } from './modules/local/trim.nf'
include { ASSEMBLE } from './modules/local/assembly.nf'
include { MAP } from './modules/local/mapping.nf'
include { COVERAGE } from './modules/local/coverage.nf'
include { INDEX_REFERENCE } from './modules/local/index.nf'
include { DOWNLOAD } from './modules/local/download.nf'
include { FILTER_VARIANTS } from './modules/local/filter_variants.nf'

include { FASTQC as FASTQC_RAW } from './modules/local/fastqc.nf'
include { FASTQC as FASTQC_TRIM } from './modules/local/fastqc.nf'

include { BCFTOOLS_MPILEUP } from './modules/nf-core/bcftools/mpileup/main.nf'

workflow CONTROL_WORKFLOW {

    take:
    reads_ch
    ref_ch
    ref_indexed

    main:

    control_qc_raw = FASTQC_RAW(reads_ch)

    control_trim = TRIM(reads_ch)

    control_qc_trim = FASTQC_TRIM(control_trim)

    control_bam = MAP(control_trim, ref_ch)

    control_cov =
        COVERAGE(
            control_bam.map { meta,bam,bai -> tuple(meta,bam) }
        )

    mpileup_ref =
        ref_indexed.map { fasta, fai ->
            tuple([id:'reference'], fasta, fai)
        }

    control_mpileup_bam =
        control_bam.map { meta, bam, bai ->
            tuple(meta, bam, [], [])
        }

    control_vcf =
        BCFTOOLS_MPILEUP(
            control_mpileup_bam,
            mpileup_ref,
            false
        )

    emit:
    vcf = control_vcf.vcf
}

workflow TREATED_WORKFLOW {

    take:
    reads_ch
    ref_ch
    ref_indexed

    main:

    treated_qc_raw = FASTQC_RAW(reads_ch)

    treated_trim = TRIM(reads_ch)

    treated_qc_trim = FASTQC_TRIM(treated_trim)

    treated_bam = MAP(treated_trim, ref_ch)

    treated_cov =
        COVERAGE(
            treated_bam.map { meta,bam,bai -> tuple(meta,bam) }
        )

    mpileup_ref =
        ref_indexed.map { fasta, fai ->
            tuple([id:'reference'], fasta, fai)
        }

    treated_mpileup_bam =
        treated_bam.map { meta,bam,bai ->
            tuple(meta, bam, [], [])
        }

    treated_vcf =
        BCFTOOLS_MPILEUP(
            treated_mpileup_bam,
            mpileup_ref,
            false
        )

    emit:
    vcf = treated_vcf.vcf
}

workflow {

    samples_ch =
        Channel
            .fromPath(params.samplesheet, checkIfExists:true)
            .splitCsv(header:true)
            .map { row ->
                tuple(
                    [
                        id: row.sample_id,
                        group: row.group
                    ],
                    row.sra
                )
            }

    // Download reads

    reads_ch = DOWNLOAD(samples_ch)

    // Split reads by group

    control_reads =
        reads_ch.filter { meta, reads ->
            meta.group == 'control'
        }

    treated_reads =
        reads_ch.filter { meta, reads ->
            meta.group == 'treated'
        }

    // Reference

    ref_ch = Channel.fromPath(params.reference)

    ref_indexed = INDEX_REFERENCE(ref_ch)

    // Run subworkflows

    control_result =
        CONTROL_WORKFLOW(
            control_reads,
            ref_ch,
            ref_indexed
        )

    treated_result =
        TREATED_WORKFLOW(
            treated_reads,
            ref_ch,
            ref_indexed
        )

    // Merge separated workflows

    all_variants =
        control_result.vcf.mix(
            treated_result.vcf
        )

    FILTER_VARIANTS(all_variants)
}

