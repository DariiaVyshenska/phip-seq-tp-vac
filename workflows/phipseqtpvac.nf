include { CUTADAPT } from '../modules/cutadapt/main'
include { KALLISTO_QUANT } from '../modules/kallisto/quant/main' 
include { PARSE_KALLISTO_OUTPUT } from '../modules/parse_kallisto_output/main'

workflow PHIPSEQTPVAC {
    take:
    samplesheet

    main:

    CUTADAPT(samplesheet)

    KALLISTO_QUANT(
        CUTADAPT.out.reads,
        Channel.value([null, file(params.kal_index_ref)]),
        params.kal_gtf,
        params.kal_chromosomes,
        params.kal_fragment_length,
        params.kal_fragment_length_sd
    ).results.map {_meta, kallisto_out_path -> [kallisto_out_path]}.collect()
    .set { collected_results_ch }

    PARSE_KALLISTO_OUTPUT(
        collected_results_ch,
        channel.fromPath(params.target_keys)
    )
}