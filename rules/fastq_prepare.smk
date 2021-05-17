import os.path

####################################
# RAW FASTQ PREPARATION
#
rule fastq_prepare_SE:
    output: fastq = "raw_fastq/{sample_name}_SE.fastq.gz"
    log:    run = "logs/{sample_name}/fastq_prepare/fastq_prepare_SE.log"
    params: in_filename = lambda wildcards: cfg.loc[cfg[SAMPLE] == wildcards.sample,'fastq_filename'].min(),
            umi = lambda wildcards: cfg.loc[cfg[SAMPLE] == wildcards.sample,'umi'].min(),
            run_name = lambda wildcards: cfg.loc[cfg[SAMPLE] == wildcards.sample,'run_name'].min(),
            lib_name = lambda wildcards: cfg.loc[cfg[SAMPLE] == wildcards.sample,'lib_name'].min(),
            is_first = lambda wildcards: cfg[SAMPLE].sort_values().tolist()[0] == wildcards.sample,
            rep_fastq = lambda wildcards: set(expand(DIR + "/raw_fastq/{sample}_SE.fastq.gz",sample = cfg.loc[cfg[SAMPLE] == wildcards.sample, SAMPLE_REP]))
    threads:  1
    conda:  "../wrappers/fastq_prepare_SE/env.yaml"
    script:  "../wrappers/fastq_prepare_SE/script.py"

rule fastq_prepare_PE:
    output:
        R1 = "raw_fastq/{sample_name}_R1.fastq.gz",
        R2 = "raw_fastq/{sample_name}_R2.fastq.gz"
    log:
        run = "logs/{sample_name}/fastq_prepare/fastq_prepare_PE.log"
    params:
        in_filename = lambda wildcards: cfg.loc[(cfg[SAMPLE] == wildcards.sample),'fastq_filename'].min(),
        umi = config["libraries"]["UMI"],
        run_name = lambda wildcards: cfg.loc[(cfg[SAMPLE] == wildcards.sample),'run_name'].min(),
        rep_fastq_R1 = "raw_fastq/{sample_name}_R1.fastq.gz",
        rep_fastq_R2 = "raw_fastq/{sample_name}_R2.fastq.gz
    threads:  1
    conda:  "../wrappers/fastq_prepare_PE/env.yaml"
    script:  "../wrappers/fastq_prepare_PE/script.py"
