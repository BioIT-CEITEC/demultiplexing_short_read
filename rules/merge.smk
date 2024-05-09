rule merge:
    output: res_file = library_output + "/raw_fastq/{filename}",
    shell:
        "cat */raw_fastq/{wildcards.filename} > {output.res_file}"


rule create_read_count_stats:
    input: fastq_files = resulting_fastq_files,
    output: library_output + "/sequencing_run_info/samplesNumberReads.json"
    params: config = config,
            lib_name = library_output,
    script: "../wrappers/create_read_count_stats/script.py"
