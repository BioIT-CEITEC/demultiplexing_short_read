rule fastq_prepare_SE:
    input:  in_filename = lambda wildcards: expand("{run_name}/{sample}_S{sample_index}_R1_001.fastq.gz",run_name = config["run_name"]\
                    ,sample = wildcards.sample\
                    ,sample_index = sample_tab.loc[(sample_tab["sample_name"] == wildcards.sample) & (sample_tab.library == wildcards.library)].index[0])[0]
    output: fastq = "{library}/raw_fastq/{sample}.fastq.gz"
    log:    "{library}/sample_logs/{sample}/fastq_prepare_SE.log"
    params: umi = lambda wildcards: config["libraries"][wildcards.library]["UMI"],
            run_name = config["run_name"],
    threads:  1
    conda:  "wrappers/fastq_prepare_SE/env.yaml"
    script: "wrappers/fastq_prepare_SE/script.py"

rule fastq_prepare_PE:
    input:  in_filename = lambda wildcards: expand("{run_name}/{sample}_S{sample_index}_R1_001.fastq.gz",run_name = config["run_name"]\
                    ,sample = wildcards.sample\
                    ,sample_index = sample_tab.loc[(sample_tab["sample_name"] == wildcards.sample) & (sample_tab.library == wildcards.library)].index[0])[0]
    output: R1 = "{library}/raw_fastq/{sample}_R1.fastq.gz",
            R2 = "{library}/raw_fastq/{sample}_R2.fastq.gz"
    log:    "{library}/sample_logs/{sample}/fastq_prepare_PE.log"
    params: umi = lambda wildcards: config["libraries"][wildcards.library]["UMI"],
            run_name = config["run_name"],
    threads:  1
    conda:  "wrappers/fastq_prepare_PE/env.yaml"
    script: "wrappers/fastq_prepare_PE/script.py"