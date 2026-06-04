nextflow.enable.dsl=2

include { TRIM } from './modules/local/trim.nf'
include { ASSEMBLE } from './modules/local/assembly.nf'
include { MAP } from './modules/local/mapping.nf'
include { COVERAGE } from './modules/local/coverage.nf'
include { INDEX_REFERENCE } from './modules/local/index.nf'
include { DOWNLOAD } from './modules/local/download.nf'
include { FASTQC as FASTQC_RAW }  from './modules/local/fastqc.nf'
include { FASTQC as FASTQC_TRIM } from './modules/local/fastqc.nf'
include { BCFTOOLS_MPILEUP } from './modules/nf-core/bcftools/mpileup/main.nf'

workflow {

    def reads_ch

    if (params.reads) {

        if (params.paired) {
            reads_ch = Channel
                .fromFilePairs(params.reads, checkIfExists: true)
                .map { id, reads -> [id, reads] }

        } else {
            reads_ch = Channel
                .fromPath(params.reads, checkIfExists: true)
                .map { file -> [file.baseName, [file]] }
        }

    } else if (params.sra_ids) {

        reads_ch = DOWNLOAD(params.sra_ids)

    } else {
        error "Provide either --reads or --sra_ids"
    }

    // QC raw reads
    qc_raw = FASTQC_RAW(reads_ch)

    // Trimming
    trimmed = TRIM(reads_ch)

    // QC trimmed
    qc_trim = FASTQC_TRIM(trimmed)

    def ref_ch

    if (params.reference) {
        ref_ch = Channel.fromPath(params.reference)
    } else {
        ref_ch = ASSEMBLE(trimmed)
    }

    // Mapping
    bam = MAP(trimmed, ref_ch)

    // Coverage
    COVERAGE(bam, ref_ch)

    // Reference indexing
    ref_indexed = INDEX_REFERENCE(ref_ch)

    mpileup_bam = bam.map { id, bamfile, bai ->
        tuple(
            [id:id],
            bamfile,
            [],
            []
        )
    }

    mpileup_ref = ref_indexed.map { fasta, fai ->
        tuple(
            [id:'reference'],
            fasta,
            fai
        )
    }

    // Variant calling using nf-core module
    vcf = BCFTOOLS_MPILEUP(
        mpileup_bam,
        mpileup_ref,
        false
    )
}