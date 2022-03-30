
rule merge:
    output: res_file = library_output + "/raw_fastq/{filename}",
    shell:
        "cat */raw_fastq/{wildcards.filename} > {output.res_file}"

rule creat_read_count_stats:
    input: all_sample_inputs
    output: library_output + "/sequencing_run_info/Stats.json"
    params: sample_tab = sample_tab
    shell:
        "cat */raw_fastq/{wildcards.filename} > {output.res_file}"
