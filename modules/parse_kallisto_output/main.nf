process PARSE_KALLISTO_OUTPUT {
    tag "kallisto_out all"
    label 'process_low'

    container 'docker.io/dariiavyshenska/python-phipseq'

    input:
    path kall_output_files
    path keys_csv

    output:
    file "kallisto_raw_counts_merged.csv"
    file "kallisto_parser.log"


    script:
    """
    kallisto_output_parser.py . $keys_csv $kall_output_files >> kallisto_parser.log 2>&1
    """
}
