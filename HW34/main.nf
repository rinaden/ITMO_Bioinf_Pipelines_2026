nextflow.enable.dsl=2

include { TRIM } from './modules/trim.nf'
include { ASSEMBLE } from './modules/assembly.nf'
include { MAP } from './modules/mapping.nf'
include { COVERAGE } from './modules/coverage.nf'
include { DOWNLOAD } from './modules/download.nf'
include { FASTQC as FASTQC_RAW }  from './modules/fastqc.nf'
include { FASTQC as FASTQC_TRIM } from './modules/fastqc.nf'
include { CALLING } from './modules/calling.nf'

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

    // Variant calling
    vcf = CALLING(bam, ref_ch)
}