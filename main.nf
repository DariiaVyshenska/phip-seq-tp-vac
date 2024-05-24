#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    greninger-lab/phipseqtpvac
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/greninger-lab/phipseqtpvac
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
 include { CUTADAPT } from './modules/nf-core/cutadapt/main'
 include { KALLISTO_QUANT } from './modules/nf-core/kallisto/quant/main' 
 include { PARSE_KALLISTO_OUTPUT } from './subworkflows/local/parse_kallisto_output/main'
// include { PHIPSEQTPVAC  } from './workflows/phipseqtpvac'
// include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_phipseqtpvac_pipeline'
// include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_phipseqtpvac_pipeline'

// include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_phipseqtpvac_pipeline'

// processes - will need to extract them in separate files:



/*

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
// workflow NFCORE_PHIPSEQTPVAC {

//     take:
//     samplesheet // channel: samplesheet read in from --input

//     main:

//     //
//     // WORKFLOW: Run pipeline
//     //
//     PHIPSEQTPVAC (
//         samplesheet
//     )

//     emit:
//     multiqc_report = PHIPSEQTPVAC.out.multiqc_report // channel: /path/to/multiqc_report.html

// }
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.indir = file(params.indir).toAbsolutePath()
params.outdir = file(params.outdir).toAbsolutePath()
params.target_keys = file(params.target_keys).toAbsolutePath()

workflow {
    // set default ref to kalisto (aws link?) & set it as optional parameter from user

    Channel.fromFilePairs("${params.indir}/*_R{1,2}.fastq.gz")
    .map { id, reads -> tuple([id: id, single_end: reads.size() == 1], reads) }
    .set{ fastq_pairs_ch }

    Channel.fromPath(params.target_keys)
    .set{ library_target_keys_csv }


    CUTADAPT(fastq_pairs_ch)

    KALLISTO_QUANT(
        CUTADAPT.out.reads,
        Channel.value([null, file(params.kal_index_ref)]),
        params.kal_gtf,
        params.kal_chromosomes,
        params.kal_fragment_length,
        params.kal_fragment_length_sd
    ).results.map {meta, kallisto_out_path -> [kallisto_out_path]}.collect()
    .map { files -> tuple([id: 'all_abundance_files'], files)}
    .set { collected_results_ch }


    PARSE_KALLISTO_OUTPUT(
       collected_results_ch,
       library_target_keys_csv
    )

    // for future dev I can incorporate generating the index too:
        // #make an index
        // kallisto index yourindex.fasta -i output_filename

    // main:

//     //
//     // SUBWORKFLOW: Run initialisation tasks
//     //
//     PIPELINE_INITIALISATION (
//         params.version,
//         params.help,
//         params.validate_params,
//         params.monochrome_logs,
//         args,
//         params.outdir,
//         params.input
//     )

//     //
//     // WORKFLOW: Run main workflow
//     //
//     NFCORE_PHIPSEQTPVAC (
//         PIPELINE_INITIALISATION.out.samplesheet
//     )

//     //
//     // SUBWORKFLOW: Run completion tasks
//     //
//     PIPELINE_COMPLETION (
//         params.email,
//         params.email_on_fail,
//         params.plaintext_email,
//         params.outdir,
//         params.monochrome_logs,
//         params.hook_url,
//         NFCORE_PHIPSEQTPVAC.out.multiqc_report
//     )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
