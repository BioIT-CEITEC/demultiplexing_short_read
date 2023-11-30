if config["run_sequencer_type"] != "AVITI":
    rule illumina_create_samplesheet:
        input:  run_info = expand("{run_dir}/RunInfo.xml",run_dir=config["run_dir"]),
        output: samplesheet_csv = config["run_name"] + "/{bcl2fastq_params_slug}/run_samplesheet.csv",
        params: sample_tab = lambda wildcards: sample_tab[sample_tab["bcl2fastq_params_slug"] == wildcards.bcl2fastq_params_slug],
                date = config["run_date"],
                run_name = config["run_name"],
                run_forward_read_length = config["run_forward_read_length"],
                run_reverse_read_length = config["run_reverse_read_length"]
        script: "../wrappers/create_samplesheet/script.py"

    rule illumina_bcl2fastq:
        input:  run_complete_check = expand("{run_dir}/RTAComplete.txt",run_dir = config["run_dir"]),
                samplesheet_csv = config["run_name"] + "/{bcl2fastq_params_slug}/run_samplesheet.csv"
        output: demultiplex_complete_check = config["run_name"] + "/{bcl2fastq_params_slug}/Reports/html/index.html",
                stats = config["run_name"] + "/{bcl2fastq_params_slug}/Stats/Stats.json",
                html  = config["run_name"] + "/{bcl2fastq_params_slug}/Stats/bcl2fastq_multiqc.html",
                mzip  = config["run_name"] + "/{bcl2fastq_params_slug}/Stats/bcl2fastq_multiqc_data.zip",
        params: library_configs = lambda wildcards: {lib_name:config["libraries"][lib_name] for lib_name in set(sample_tab[sample_tab["bcl2fastq_params_slug"] == wildcards.bcl2fastq_params_slug].library)},
                tmp_dir = GLOBAL_TMPD_PATH
        log:    config["run_name"] + "/{bcl2fastq_params_slug}/bcl2fastq.log"
        params: library_configs = lambda wildcards: {lib_name:config["libraries"][lib_name] for lib_name in set(sample_tab[sample_tab["bcl2fastq_params_slug"] == wildcards.bcl2fastq_params_slug].library)}
        conda: "../wrappers/bcl2fastq/env.yaml"
        script: "../wrappers/bcl2fastq/script.py"

    rule illumina_stats_copy:
        input:  stats = lambda wildcards: expand(config["run_name"] + "/{bcl2fastq_params_slug}/Stats/Stats.json", bcl2fastq_params_slug = sample_tab.loc[sample_tab.library == wildcards.library_name,'bcl2fastq_params_slug'].min())[0],
                html  = lambda wildcards: expand(config["run_name"] + "/{bcl2fastq_params_slug}/Stats/bcl2fastq_multiqc.html", bcl2fastq_params_slug = sample_tab.loc[sample_tab.library == wildcards.library_name,'bcl2fastq_params_slug'].min())[0],
                mzip  = lambda wildcards: expand(config["run_name"] + "/{bcl2fastq_params_slug}/Stats/bcl2fastq_multiqc_data.zip", bcl2fastq_params_slug = sample_tab.loc[sample_tab.library == wildcards.library_name,'bcl2fastq_params_slug'].min())[0],
        output: stats = "{library_name}/sequencing_run_info/Stats.json",
                html  = "{library_name}/sequencing_run_info/bcl2fastq_multiqc.html",
                mzip  = "{library_name}/sequencing_run_info/bcl2fastq_multiqc_data.zip",
        log:    run = "{library_name}/sequencing_run_info/stats_copy.log",
        script: "../wrappers/copy_stats/script.py"

    rule illumina_fastq_mv:
        input:  demultiplex_complete_check = expand(config["run_name"] + "/{bcl2fastq_params_slug}/Reports/html/index.html",bcl2fastq_params_slug = list(pd.unique(sample_tab['bcl2fastq_params_slug']))),
        output: fastqs_out = expand(config["run_name"] + "/{sample}_S{sample_index}_R1_001.fastq.gz",zip,sample = sample_tab.sample_name\
                                                                                ,sample_index = range(1,len(sample_tab.index)+1))
        params: fastqs_in = expand(config["run_name"] + "/{bcl2fastq_params_slug}/{sample}_S{sample_index}_R1_001.fastq.gz",zip\
                                                                                ,sample = sample_tab.sample_name\
                                                                                ,sample_index = sample_tab.slug_id\
                                                                                ,bcl2fastq_params_slug = sample_tab.bcl2fastq_params_slug)
        script: "../wrappers/fastq_mv/script.py"

    rule illumina_fastq_mv_to_lib:
        input:  lambda wildcards: expand("{run_name}/{sample_name}_S{sample_index}_{read_num}_001.fastq.gz",run_name = config["run_name"]\
                        ,sample_name = wildcards.sample_name \
                        ,sample_index = sample_tab.loc[(sample_tab["sample_name"] == wildcards.sample) & (sample_tab.library == wildcards.library)].index[0]\
                        ,read_num = wildcards.read_num)[0]
        output: "{library}/raw_fastq/{sample_name}_{read_num}.fastq.gz",
        shell:
            "mv {input} {output}"