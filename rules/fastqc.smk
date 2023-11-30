

#TODO finish adaptor check running must be per library only for libraries with the checkd parameter
# def minion_input(wcs):
#     if config['check_adaptors']:
#         "qc_reports/raw_fastq_minion_adaptors_mqc.tsv"

rule merge_fastq_qc:
    input:  html = all_sample_inputs
    output: html = "qc_reports/raw_fastq_multiqc.html"
    log:    "logs/merge_fastq_qc.log"
    conda:  "../wrappers/merge_fastq_qc/env.yaml"
    script: "../wrappers/merge_fastq_qc/script.py"

rule raw_fastq_qc:
    input:  raw_fastq = "{library}/raw_fastq/{sample_name}_{read_num}.fastq.gz",
    output: html = "{library}/qc_reports/{sample_name}/raw_fastqc/{sample_name}_{read_num}_fastqc.html"
    log:    "{library}/logs/{sample_name}/raw_fastqc_{read_num}.log"
    params: extra = "--noextract --format fastq --nogroup",
    threads:  1
    conda:  "../wrappers/raw_fastq_qc/env.yaml"
    script: "../wrappers/raw_fastq_qc/script.py"
