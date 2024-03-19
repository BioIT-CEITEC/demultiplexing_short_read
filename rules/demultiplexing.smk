if config["run_sequencer_type"] == "AVITI":
    rule aviti_create_samplesheet:
        input:  run_info = expand("{run_dir}/RunParameters.json",run_dir=config["run_dir"]),
        output: run_manifest = "run_manifest.csv",
        params: sample_tab = sample_tab,
                config = config
        script: "../wrappers/aviti_create_samplesheet/script.py"

    rule aviti_Bases2Fastq:
        input:  run_complete_check = expand("{run_dir}/RunUploaded.json",run_dir = config["run_dir"]),
                run_manifest = "run_manifest.csv"
        output: fastq_files = expand("Samples/{library}/{sample_name}/{sample_name}_R{read_num}.fastq.gz",zip,library = sample_file_tab.library\
                                                                                                      ,sample_name=sample_file_tab.sample_name\
                                                                                                      ,read_num=sample_file_tab.read_num),
                stats = expand("Samples/{library}/{library}_{run_info_suffix}",library = sample_file_tab.library
                                                                            ,run_info_suffix = ["QC.html","Metrics.csv","RunStats.json","IndexAssignment.csv"])
                # demultiplex_complete_check = config["run_name"] + "/{bcl2fastq_params_slug}/Reports/html/index.html",
                # stats = config["run_name"] + "/{bcl2fastq_params_slug}/Stats/Stats.json",
                # html  = config["run_name"] + "/{bcl2fastq_params_slug}/Stats/bcl2fastq_multiqc.html",
                # mzip  = config["run_name"] + "/{bcl2fastq_params_slug}/Stats/bcl2fastq_multiqc_data.zip",
        params: tmp_dir = GLOBAL_TMPD_PATH
        threads: 30
        log:    "logs/Bases2Fastq.log"
        params: library_configs = lambda wildcards: {lib_name:config["libraries"][lib_name] for lib_name in set(sample_tab[sample_tab["bcl2fastq_params_slug"] == wildcards.bcl2fastq_params_slug].library)}
        conda: "../wrappers/aviti_bases2fastq/env.yaml"
        script: "../wrappers/aviti_bases2fastq/script.py"

    rule aviti_fastq_mv:
        input:  "Samples/{library}/{sample_name}/{sample_name}_R{read_num}.fastq.gz"
        output: "{library}/raw_fastq/{sample_name}_R{read_num}.fastq.gz",
        shell:
            "mkdir -p '$(dirname {output})'; mv {input} {output}"

    rule aviti_stats_copy:
        input:  files = expand("Samples/{{library}}/{{library}}_{run_info_suffix}",run_info_suffix = ["QC.html","Metrics.csv","RunStats.json","IndexAssignment.csv"])
        output: stats = "{library}/sequencing_run_info/Stats.json",
        log:    run = "{library}/sequencing_run_info/stats_copy.log",
        params: stats_json_file = "Samples/{library}/{library}_RunStats.json"
        script: "../wrappers/aviti_stats_copy/script.py"

elif config["run_sequencer_type"] == "MGI":
    rule mgi_create_samplesheet:
        input: run_info=expand("{run_dir}/L01/BioInfo.csv",run_dir=config["run_dir"]),
        output: run_manifest="mgi_sample_sheet_{lane}.txt"
        params: sample_tab=sample_tab,
            config=config
        script: "../wrappers/mgi_create_samplesheet/script.py"

    rule mgi_calDemux:
        input:  run_manifest="mgi_sample_sheet_{lane}.txt"
        output: fastq_files = expand("demux/FS2000/{{lane}}/FS2000_{{lane}}_{sample_name}_{read_num}.fq.gz",zip \
                                            ,sample_name=sample_file_tab.sample_name \
                                            ,read_num=sample_file_tab.read_num),
                stats=expand("demux/FS2000/{{lane}}/{stat_filenames}.txt" \
                                            ,stat_filenames=["SequenceStat", "BarcodeStat"])
        # demultiplex_complete_check = config["run_name"] + "/{bcl2fastq_params_slug}/Reports/html/index.html",
        # stats = config["run_name"] + "/{bcl2fastq_params_slug}/Stats/Stats.json",
        # html  = config["run_name"] + "/{bcl2fastq_params_slug}/Stats/bcl2fastq_multiqc.html",
        # mzip  = config["run_name"] + "/{bcl2fastq_params_slug}/Stats/bcl2fastq_multiqc_data.zip",
        params: tmp_dir=GLOBAL_TMPD_PATH
        threads: 30
        log: "logs/calDemux_{lane}.log"
        conda: "../wrappers/mgi_calDemux/env.yaml"
        script: "../wrappers/mgi_calDemux/script.py"


    def mgi_cat_lane_fastq_input(wildcards):
        return expand("demux/FS2000/{lane}/FS2000_{lane}_"+wildcards.sample_name+"_"+wildcards.read_num+".fq.gz" \
            ,lane = per_library_used_lanes[wildcards.library])

    rule mgi_cat_lane_fastq:
        input: mgi_cat_lane_fastq_input
        output: "{library}/raw_fastq/{sample_name}_R{read_num}.fastq.gz",
        threads: 30
        shell:
            "mkdir -p '$(dirname {output})'; cat {input} > {output}"


    def mgi_stats_copy_input(wildcards):
        return expand("demux/FS2000/{lane}/{stat_filenames}.txt" \
            ,lane = per_library_used_lanes[wildcards.library]
            ,stat_filenames = ["SequenceStat", "BarcodeStat"])

    rule mgi_stats_copy:
        input: mgi_stats_copy_input
        output: "{library}/sequencing_run_info/Stats.json",
        log: run="{library}/sequencing_run_info/stats_copy.log",
        params: stats_json_file="Samples/{library}/{library}_RunStats.json"
        # script: "../wrappers/mgi_stats_copy/script.py"
        shell:
            "mkdir -p '$(dirname {output})'; cat {input} > {output}"

    # rule illumina_fastq_mv:
    #     input:  demultiplex_complete_check = expand(config["run_name"] + "/{bcl2fastq_params_slug}/Reports/html/index.html",bcl2fastq_params_slug = list(pd.unique(sample_tab['bcl2fastq_params_slug']))),
    #     output: fastqs_out = expand(config["run_name"] + "/{sample}_S{sample_index}_R1_001.fastq.gz",zip,sample = sample_tab.sample_name\
    #                                                                             ,sample_index = range(1,len(sample_tab.index)+1))
    #     params: fastqs_in = expand(config["run_name"] + "/{bcl2fastq_params_slug}/{sample}_S{sample_index}_R1_001.fastq.gz",zip\
    #                                                                             ,sample = sample_tab.sample_name\
    #                                                                             ,sample_index = sample_tab.slug_id\
    #                                                                             ,bcl2fastq_params_slug = sample_tab.bcl2fastq_params_slug)
    #     script: "../wrappers/fastq_mv/script.py"
    #
    # rule illumina_fastq_mv_to_lib:
    #     input:  lambda wildcards: expand("{run_name}/{sample_name}_S{sample_index}_{read_num}_001.fastq.gz",run_name = config["run_name"]\
    #                     ,sample_name = wildcards.sample_name \
    #                     ,sample_index = sample_tab.loc[(sample_tab["sample_name"] == wildcards.sample) & (sample_tab.library == wildcards.library)].index[0]\
    #                     ,read_num = wildcards.read_num)[0]
    #     output: "{library}/raw_fastq/{sample_name}_{read_num}.fastq.gz",
    #     shell:
    #         "mv {input} {output}"
else:
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