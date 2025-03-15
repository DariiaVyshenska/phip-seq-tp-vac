include { validateParameters } from 'plugin/nf-validation'

workflow UTILS_NFVALIDATION_PLUGIN {

    take:
    validate_params  // boolean: validate parameters
    schema_filename  //    path: JSON schema file, null to use default value

    main:

    log.debug "Using schema file: ${schema_filename}"

    if (validate_params){
        validateParameters(parameters_schema: schema_filename)
    }

    emit:
    dummy_emit = true
}
