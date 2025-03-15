//
// Subworkflow with utility functions specific to the nf-core pipeline template
//

import org.yaml.snakeyaml.Yaml // import org.yaml.snakeyaml.Yaml: This line imports the SnakeYAML library, which is used for parsing and generating YAML files in Java/Groovy. In the context of nf-core pipelines, it's likely used for reading and processing YAML configuration files or for generating YAML-formatted output (e.g., for version information).
import nextflow.extension.FilesEx // import nextflow.extension.FilesEx: This imports Nextflow's extended file handling utilities. FilesEx provides additional file-related operations that may not be available in standard Java file handling classes.

/*
========================================================================================
    SUBWORKFLOW DEFINITION
========================================================================================
*/

workflow UTILS_NFCORE_PIPELINE {

    take:
    nextflow_cli_args

    main:
    valid_config = checkConfigProvided()
    checkProfileProvided(nextflow_cli_args)

    emit:
    valid_config
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
//  Warn if a -profile or Nextflow config has not been provided to run the pipeline
//
def checkConfigProvided() {
    valid_config = true
    if (workflow.profile == 'standard' && workflow.configFiles.size() <= 1) {
        log.warn "[$workflow.manifest.name] You are attempting to run the pipeline without any custom configuration!\n\n" +
            "This will be dependent on your local compute environment but can be achieved via one or more of the following:\n" +
            "   (1) Using an existing pipeline profile e.g. `-profile docker` or `-profile singularity`\n" +
            "   (2) Using an existing nf-core/configs for your Institution e.g. `-profile crick` or `-profile uppmax`\n" +
            "   (3) Using your own local custom config e.g. `-c /path/to/your/custom.config`\n\n" +
            "Please refer to the quick start section and usage docs for the pipeline.\n "
        valid_config = false
    }
    return valid_config
}

//
// Exit pipeline if --profile contains spaces
//
def checkProfileProvided(nextflow_cli_args) {
    if (workflow.profile.endsWith(',')) {
        error "The `-profile` option cannot end with a trailing comma, please remove it and re-run the pipeline!\n" +
            "HINT: A common mistake is to provide multiple values separated by spaces e.g. `-profile test, docker`.\n"
    }
    if (nextflow_cli_args[0]) {
        log.warn "nf-core pipelines do not accept positional arguments. The positional argument `${nextflow_cli_args[0]}` has been detected.\n" +
            "HINT: A common mistake is to provide multiple values separated by spaces e.g. `-profile test, docker`.\n"
    }
}


// can I reulse it?? future-dev : concat workflow vesions into a file
// two imports at the top are probably for these functions.
//
// Generate workflow version string
//
def getWorkflowVersion() {
    String version_string = ""
    if (workflow.manifest.version) {
        def prefix_v = workflow.manifest.version[0] != 'v' ? 'v' : ''
        version_string += "${prefix_v}${workflow.manifest.version}"
    }

    if (workflow.commitId) {
        def git_shortsha = workflow.commitId.substring(0, 7)
        version_string += "-g${git_shortsha}"
    }

    return version_string
}

//
// Get software versions for pipeline
//
def processVersionsFromYAML(yaml_file) {
    Yaml yaml = new Yaml()
    versions = yaml.load(yaml_file).collectEntries { k, v -> [ k.tokenize(':')[-1], v ] }
    return yaml.dumpAsMap(versions).trim()
}

//
// Get workflow version for pipeline
//
def workflowVersionToYAML() {
    return """
    Workflow:
        $workflow.manifest.name: ${getWorkflowVersion()}
        Nextflow: $workflow.nextflow.version
    """.stripIndent().trim()
}

//
// Get channel of software versions used in pipeline in YAML format
//
def softwareVersionsToYAML(ch_versions) {
    return ch_versions
                .unique()
                .map { processVersionsFromYAML(it) }
                .unique()
                .mix(Channel.of(workflowVersionToYAML()))
}
