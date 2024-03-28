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
// include { PHIPSEQTPVAC  } from './workflows/phipseqtpvac'
// include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_phipseqtpvac_pipeline'
// include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_phipseqtpvac_pipeline'

// include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_phipseqtpvac_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
// params.fasta = getGenomeAttribute('fasta')

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
// params.str = "/Users/mac_studio_np/Documents/nextflow_repos/phipseq/nf-core-phipseqtpvac"

// process printInputFilenames {
//     input:
//     path fastaFilesDir

//     output:
//     stdout

//     """
//     printf "%s\n" "Hello from the main workflow"
//     // printf '${params.str}'
//     """
// }


// need to make sure both absolute and relative paths are taken into account:
// path types: ./data, data, /my/abs/path/data

// process showFilesInDir {
//     // publishDir '/Users/mac_studio_np/Documents/nextflow_repos/phipseq/nf-core-phipseqtpvac/outputDir', mode: 'copy'

//     // output:
//     // path "files_ls_out.txt"
//     output:
//     stdout

//     script:
//     def inputFastaDirABS = new File(params.inputFastaDir).getAbsolutePath()
//     // """
//     // ls -la '${inputFastaDirABS}' > files_ls_out.txt
//     // """
//     """
//     ls -la '${inputFastaDirABS}'
//     """
// }


// val a = Channel.value([null, file(params.kal_index_ref)])
workflow {
    // make sure to grab unzipped files if it's appropriate!
    def inputFastqDirABS = new File(params.inputFastaDir).getAbsolutePath()


    Channel.fromFilePairs("${inputFastqDirABS}/**/*_L001_R{1,2}_001.fastq.gz")
    .map { id, reads -> tuple([id: id, single_end: reads.size() == 1], reads) }
    .set{ fastq_pairs_ch }

// RUNNING FROM COMMAND LINE
// nextflow run main.nf --inputFastaDir raw_test_data -profile docker -c ~/nextflow_configs/learning.config; mv .nextflow.log* ./outputDir/pipeline_inf

    CUTADAPT(fastq_pairs_ch)


    // kallisto quant 
    //     --index ./ref/tp_oligos_kallisto             => index file
    //     -o ./tp_abundance/${sample}_kallisto         => output dir
    //     --plaintext                                  => Output plaintext instead of HDF5
    //     ./trim/${sample}_R1.fastq.gz                 => ...reads
    //     ./trim/${sample}_R2.fastq.gz 

// "kallisto quant 
// --threads 10 
// --index tp_oligos_kallisto 
// --plaintext 
// -o 240215_beads_1_S5 
// 240215_beads_1_S5_1.trim.fastq.gz 
// 240215_beads_1_S5_2.trim.fastq.gz"

    // mkdir -p $prefix && kallisto quant \\
    //     --threads ${task.cpus} \\
    //     --index ${index} \\
    //     ${gtf_input} \\
    //     ${chromosomes_input} \\
    //     ${single_end_params} \\
    //     ${strandedness} \\
    //     ${args} \\
    //     -o $prefix \\
    //     ${reads} 2> >(tee -a ${prefix}/kallisto_quant.log >&2)
    
    KALLISTO_QUANT(
        CUTADAPT.out.reads,
        Channel.value([null, file(params.kal_index_ref)]),
        "",
        ""
        // false,
        // false,
        // "",
        // ""
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
