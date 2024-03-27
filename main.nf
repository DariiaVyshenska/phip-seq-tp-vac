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
params.inputFastaDir = './data'

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



workflow {
    // showFilesInDir | view { it }

    // make sure to grab unzipped files if it's appropriate!
    def inputFastqDirABS = new File(params.inputFastaDir).getAbsolutePath()
    // fastqFileChannel = Channel.fromPath("${inputFastqDirABS}/**/*.fastq.gz")

    // fastqFileChannel = Channel.fromPath("${inputFastqDirABS}/**/*_L001_R{1,2}_001.fastq.gz")
    Channel.fromFilePairs("${inputFastqDirABS}/**/*_L001_R{1,2}_001.fastq.gz")
    .map { id, files -> tuple([id: id, single_end: files.size() == 1], files) }
    .set{ fastq_pairs_ch }

// RUNNING FROM COMMAND LINE
// nextflow run main.nf --inputFastaDir raw_test_data -profile docker -c ~/nextflow_configs/learning.config; mv .nextflow.log* ./outputDir/pipeline_inf


    // fastq_pairs_ch.view()
    CUTADAPT(fastq_pairs_ch)

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
