import os

## ANNOTATION of VARIANTS in SAMPLES
rule create_samplesheet:
    params: sample_tab = sample_tab
    output: samplesheet_csv = "run_samplesheet.csv",
    shell:  "touch {output.samplesheet_csv}"

rule bcl2fastq:
    input:  run_complete_check = expand("{run_dir}/RTAComplete.txt",run_dir = config["run_dir"]),
            samplesheet_csv= "run_samplesheet.csv"
    params: library_list = set(sample_tab.library)
    output: expand("{prefix}_{read_pair_tag}.fastq.gz",prefix = expand("{library}/raw_fastq/{sample}",zip,sample = sample_tab.sample_name,library = sample_tab.library),read_pair_tag = read_pair_tags),
    log:    "logs/bcl2fastq.log"
    shell:  "mkdir -p {params.library_list};touch {output}"
