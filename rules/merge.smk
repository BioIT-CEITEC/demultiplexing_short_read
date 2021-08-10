
rule merge:
    output: res_file = config["libraries"].keys()[0] + "-merged/raw_fastq/{filename}",
    shell:
        "cat */raw_fastq/{wildcards.filename} > {output.res_file}"



