params.input_reads_folder = 'data'

process run_qc {
    input:
        tuple val(reads_type), val(reads_label), path(reads)
    output:
        path "${reads_type}_qc_report/", type: 'folder'
    script:
    """
    mkdir ${reads_type}_qc_report
    fastqc -o ${reads_type}_qc_report/ $reads
    """
}

process trimm {
    input:
        tuple val(reads_label), path(reads)
    output:
        tuple val(reads_label), path("out_R?_p.fq.gz")
    script:
    """
    trimmomatic PE $reads \
        out_R1_p.fq.gz out_R1_u.fq.gz \
        out_R2_p.fq.gz out_R2_u.fq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36
    """
}

workflow {

    input_reads = channel
        .fromFilePairs("${params.input_reads_folder}/*_{1,2}.fq")

    // Initial QC input
    qc_initial_input = input_reads.map { label, files ->
        tuple('initial', label, files)
    }

    // Trimming
    trimm_result = trimm(input_reads)

    // Trimmed QC input
    qc_trimmed_input = trimm_result.map { label, files ->
        tuple('trimmed', label, files)
    }

    qc_all_input = qc_initial_input.mix(qc_trimmed_input)

    qc_results = run_qc(qc_all_input)

    publish:
        trimmed_reads = trimm_result
        qc = qc_results
}

output {
    trimmed_reads {
        path 'trimmed_reads'
    }
    qc {
        path 'qc'
    }
}