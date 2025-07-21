

#TODO finish adaptor check running must be per library only for libraries with the checkd parameter
# def minion_input(wcs):
#     if config['check_adaptors']:
#         "qc_reports/raw_fastq_minion_adaptors_mqc.tsv"

def merge_fastq_qc_input(wildcards):
    if "merged" in config and config["merged"]:
        return all_merge_sample_fastqc_files
    else:
        return expand(wildcards.library + "/qc_reports/{sample_name}/raw_fastqc/{sample_name}_R{read_num}_fastqc.html", zip \
                ,sample_name=sample_file_tab.loc[sample_file_tab.library == wildcards.library].sample_name \
                ,read_num=sample_file_tab.loc[sample_file_tab.library == wildcards.library].read_num)

rule merge_fastq_qc:
    input:  html = merge_fastq_qc_input
    output: html = "{library}/qc_reports/raw_fastq_multiqc.html"
    log:    "{library}/logs/merge_fastq_qc.log"
    conda:  "../wrappers/merge_fastq_qc/env.yaml"
    script: "../wrappers/merge_fastq_qc/script.py"

rule raw_fastq_qc:
    input:  raw_fastq = "{library}/raw_fastq/{sample_name}_R{read_num}.fastq.gz",
    output: html = "{library}/qc_reports/{sample_name}/raw_fastqc/{sample_name}_R{read_num}_fastqc.html"
    log:    "{library}/logs/{sample_name}/raw_fastqc_{read_num}.log"
    params: extra = "--noextract --format fastq --nogroup",
    threads:  1
    conda:  "../wrappers/raw_fastq_qc/env.yaml"
    script: "../wrappers/raw_fastq_qc/script.py"
