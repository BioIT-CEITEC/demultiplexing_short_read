
rule create_samplesheet:
    params: sample_tab = lambda wildcards: sample_tab[sample_tab["bcl2fastq_params_slug"] == wildcards.bcl2fastq_params_slug],
            date = config["run_date"],
            run_name = config["run_name"],
            run_forward_read_length = config["run_forward_read_length"],
            run_reverse_read_length = config["run_reverse_read_length"]
    output: samplesheet_csv = config["run_name"] + "/{bcl2fastq_params_slug}/run_samplesheet.csv",
    script: "../wrappers/create_samplesheet/script.py"


rule bcl2fastq:
    input:  run_complete_check = expand("{run_dir}/RTAComplete.txt",run_dir = config["run_dir"]),
            samplesheet_csv= config["run_name"] + "/{bcl2fastq_params_slug}/run_samplesheet.csv"
    output: demultiplex_complete_check = config["run_name"] + "/{bcl2fastq_params_slug}/Reports/html/index.html",
    params: library_configs = lambda wildcards: {lib_name:config["libraries"][lib_name] for lib_name in set(sample_tab[sample_tab["bcl2fastq_params_slug"] == wildcards.bcl2fastq_params_slug].library)}
    log:    config["run_name"] + "/{bcl2fastq_params_slug}/bcl2fastq.log"
    conda: "../wrappers/bcl2fastq/env.yaml"
    script: "../wrappers/bcl2fastq/script.py"


rule fastq_mv:
    input:  demultiplex_complete_check = expand(config["run_name"] + "/{bcl2fastq_params_slug}/Reports/html/index.html",bcl2fastq_params_slug = list(pd.unique(sample_tab['bcl2fastq_params_slug']))),
    output: fastqs_out = expand(config["run_name"] + "/{sample}_S{sample_index}_R1_001.fastq.gz",zip,sample = sample_tab.sample_name\
                                                                            ,sample_index = range(1,len(sample_tab.index)+1))
    params: fastqs_in = expand(config["run_name"] + "/{bcl2fastq_params_slug}/{sample}_S{sample_index}_R1_001.fastq.gz",zip\
                                                                            ,sample = sample_tab.sample_name\
                                                                            ,sample_index = sample_tab.slug_id\
                                                                            ,bcl2fastq_params_slug = sample_tab.bcl2fastq_params_slug)
    script: "../wrappers/fastq_mv/script.py"
