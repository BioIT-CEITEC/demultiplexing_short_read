
rule merge:
    output: res_file = list(config["library_output"].keys())[0] + "/raw_fastq/{filename}",
    shell:
        "cat */raw_fastq/{wildcards.filename} > {output.res_file}"
