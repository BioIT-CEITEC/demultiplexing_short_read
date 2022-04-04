
rule merge:
    output: res_file = library_output + "/raw_fastq/{filename}",
    shell:
        "cat */raw_fastq/{wildcards.filename} > {output.res_file}"


rule create_read_count_stats:
    input: all_sample_inputs,
    output: library_output+"/sequencing_run_info/Stats.json"
    params: config = config,
            library_output = library_output,
            lib_name = library_output+"/raw_fastq/"
    script: "../wrappers/create_read_count_stats/script.py"
