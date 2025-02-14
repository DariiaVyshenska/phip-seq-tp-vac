/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include { MULTIQC                } from '../modules/nf-core/multiqc/main'
// include { paramsSummaryMap       } from 'plugin/nf-validation'
// include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
// include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
// include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_phipseqtpvac_pipeline'

include { CUTADAPT } from '../modules/nf-core/cutadapt/main'
include { KALLISTO_QUANT } from '../modules/nf-core/kallisto/quant/main' 
include { PARSE_KALLISTO_OUTPUT } from '../subworkflows/local/parse_kallisto_output/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PHIPSEQTPVAC {
    take:
    samplesheet

    main:

    Channel.fromPath(params.target_keys)
    .set{ library_target_keys_csv }

    CUTADAPT(samplesheet)

    KALLISTO_QUANT(
        CUTADAPT.out.reads,
        Channel.value([null, file(params.kal_index_ref)]),
        params.kal_gtf,
        params.kal_chromosomes,
        params.kal_fragment_length,
        params.kal_fragment_length_sd
    ).results.map {_meta, kallisto_out_path -> [kallisto_out_path]}.collect()
    .map { files -> tuple([id: 'all_abundance_files'], files)}
    .set { collected_results_ch }


    PARSE_KALLISTO_OUTPUT(
        params.projectRoot,
        collected_results_ch,
        library_target_keys_csv
    )
}

// workflow PHIPSEQTPVAC {

//     take:
//     ch_samplesheet // channel: samplesheet read in from --input

//     main:

//     ch_versions = Channel.empty()
//     ch_multiqc_files = Channel.empty()

//     //
//     // MODULE: Run FastQC
//     //
//     FASTQC (
//         ch_samplesheet
//     )
//     ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
//     ch_versions = ch_versions.mix(FASTQC.out.versions.first())

//     //
//     // Collate and save software versions
//     //
//     softwareVersionsToYAML(ch_versions)
//         .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
//         .set { ch_collated_versions }

//     //
//     // MODULE: MultiQC
//     //
//     ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
//     ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
//     ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
//     summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
//     ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
//     ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
//     ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
//     ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
//     ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
//     ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

//     MULTIQC (
//         ch_multiqc_files.collect(),
//         ch_multiqc_config.toList(),
//         ch_multiqc_custom_config.toList(),
//         ch_multiqc_logo.toList()
//     )

//     emit:
//     multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
//     versions       = ch_versions                 // channel: [ path(versions.yml) ]
// }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
