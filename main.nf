#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    greninger-lab/phipseqtpvac
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/greninger-lab/phipseqtpvac
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2


include { PIPELINE_INITIALISATION } from './subworkflows/local/utils'
include { PHIPSEQTPVAC  } from './workflows/phipseqtpvac'

workflow {

    main:

    PIPELINE_INITIALISATION (
        params.validate_params,
        args
    )

    PHIPSEQTPVAC (
        PIPELINE_INITIALISATION.out.samplesheet
    )
}

