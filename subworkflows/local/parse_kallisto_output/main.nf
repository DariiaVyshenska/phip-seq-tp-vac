process PARSE_KALLISTO_OUTPUT {
    tag "$meta.id"
    label 'process_low'

    container 'docker.io/dariiavyshenska/python-phipseq'

    input:
    tuple val(meta), path(all_paths)
    path keys_csv

    output:
    file "kallisto_raw_counts_merged.csv"
    file "kallisto_parser.log"


    script:
    """
    python ${projectDir}/bin/kallisto_input_parser.py . ${keys_csv} ${all_paths} >> kallisto_parser.log 2>&1
    """
}
