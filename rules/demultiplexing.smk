def get_lanes_for_library(library):
    if config["run_lane_splitting"] != None:
        lane_usage = []
        for lane in range(1,config["run_lane_splitting"] + 1):
            lane_key = f"lane{lane}"
            if config["run_sequencer_type"] == "MGI":
                lane_code = f"L0{lane}"
            else:
                lane_code = f"_L00{lane}"
            if sample_tab[sample_tab['library'] == library].iloc[0][lane_key]:
                lane_usage.append(lane_code)
    else:
        if config["run_sequencer_type"] == "MGI":
            lane_usage = [entry for entry in os.listdir(config["run_dir"])
                          if os.path.isdir(os.path.join(config["run_dir"],entry)) and entry.startswith("L0")]
        else:
            lane_usage = ""
    return lane_usage


def get_demux_for_library(library):
    demux_setting = sample_tab.loc[sample_tab['library'] == library,].iloc[0]['demux_setting']
    return demux_setting


def get_sample_index_for_library(library,sample_name):
    sample_index = sample_tab.loc[(sample_tab['library'] == library) & (sample_tab['sample_name'] == sample_name)].iloc[0]['sample_index']
    return sample_index

def fastq_mv_ready_input(wildcards):
    if config["run_sequencer_type"] == "AVITI":
        input_list = expand("{demux_setting}/demux_ready.txt" \
            ,demux_setting=get_demux_for_library(wildcards.library))
    elif config["run_sequencer_type"] == "MGI":
        input_list = expand("{demux_setting}/FS2000/{lane}/demux_ready.txt",zip \
            ,lane=get_lanes_for_library(wildcards.library)
            ,demux_setting=get_demux_for_library(wildcards.library))
    else:
        input_list = expand("{demux_setting}/demux_ready.txt" \
            ,demux_setting=get_demux_for_library(wildcards.library))

    return input_list


def fastq_mv_fastq_input(wildcards):
    if config["run_sequencer_type"] == "AVITI":
        input_list = expand("{demux_setting}/Samples/"+wildcards.library+"/"+wildcards.sample_name+"/"+wildcards.sample_name+"_R"+wildcards.read_num+".fastq.gz"\
            ,demux_setting=get_demux_for_library(wildcards.library))
    elif config["run_sequencer_type"] == "MGI":
        input_list = expand("{demux_setting}/FS2000/{lane}/FS2000_{lane}_" + wildcards.sample_name + "_" + wildcards.read_num + ".fq.gz",zip \
            ,lane=get_lanes_for_library(wildcards.library)
            ,demux_setting=get_demux_for_library(wildcards.library))
    else:
        input_list = expand("{demux_setting}/" + wildcards.sample_name + "_S{sample_index}{lane}_R" + wildcards.read_num + "_001.fastq.gz",zip \
            ,lane=get_lanes_for_library(wildcards.library)
            ,demux_setting=get_demux_for_library(wildcards.library)
            ,sample_index=get_sample_index_for_library(wildcards.library,wildcards.sample_name))
        "{run_name}/{sample_name}_S{sample_index}_{read_num}_001.fastq.gz"
    return input_list

def stats_copy_input(wildcards):
    if config["run_sequencer_type"] == "AVITI":
        input_list = expand("{demux_setting}/{stat_filenames}" \
            ,demux_setting=get_demux_for_library(wildcards.library)
            ,stat_filenames=["*.html", "*.json","*.csv","*.json"])
    elif config["run_sequencer_type"] == "MGI":
        input_list = expand(expand("{demux_setting}/FS2000/{lane}/",zip \
            ,lane=get_lanes_for_library(wildcards.library)
            ,demux_setting=get_demux_for_library(wildcards.library)) + "{stat_filenames}"\
            ,stat_filenames=["SequenceStat.txt", "BarcodeStat.txt","mgi_sample_sheet*.txt"])
    else:
        input_list = expand("{demux_setting}/{stat_filenames}" \
            ,demux_setting=get_demux_for_library(wildcards.library)
            ,stat_filenames=["Reports", "Stats","run_samplesheet.csv","*.html"])

    return input_list



if config["run_sequencer_type"] == "AVITI":
    rule aviti_create_samplesheet:
        input:  run_info = expand("{run_dir}/RunParameters.json",run_dir=config["run_dir"]),
        output: run_manifest = "{demux_setting}/run_manifest.csv",
        params: config = config,
                sample_tab= lambda wildcards: sample_tab[sample_tab['demux_setting'] == wildcards.demux_setting],
        script: "../wrappers/aviti_create_samplesheet/script.py"

    rule aviti_Bases2Fastq:
        input:  run_manifest = "{demux_setting}/run_manifest.csv"
        output: demultiplex_complete = "{demux_setting}/demux_ready.txt",
                # stats = lambda wildcards: expand("{demux_setting}/Samples/{library}/{library}_{run_info_suffix}" \
                #         ,demux_setting = wildcards.demux_setting \
                #         ,library = sample_tab[sample_tab['demux_setting'] == wildcards.demux_setting,'library'].unique().tolist() \
                #         ,run_info_suffix = ["QC.html","Metrics.csv","RunStats.json","IndexAssignment.csv"])
        params: tmp_dir = GLOBAL_TMPD_PATH,
                run_dir=config["run_dir"],
                sample_tab= lambda wildcards: sample_tab[sample_tab['demux_setting'] == wildcards.demux_setting],
        threads: 30
        log:    "logs/{demux_setting}_Bases2Fastq.log"
        script: "../wrappers/aviti_bases2fastq/script.py"

    # rule aviti_stats_copy:
    #     input:  files = expand("Samples/{{library}}/{{library}}_{run_info_suffix}",run_info_suffix = ["QC.html","Metrics.csv","RunStats.json","IndexAssignment.csv"])
    #     output: stats = "{library}/sequencing_run_info/Stats.json",
    #     log:    run = "{library}/sequencing_run_info/stats_copy.log",
    #     params: stats_json_file = "Samples/{library}/{library}_RunStats.json"
    #     script: "../wrappers/aviti_stats_copy/script.py"

elif config["run_sequencer_type"] == "MGI":
    rule mgi_create_samplesheet:
        input: run_info=expand("{run_dir}/L01/BioInfo.csv",run_dir=config["run_dir"])[0],
        output: sample_sheet="mgi_sample_sheet_{lane}_{demux}.txt"
        params: lane="{lane}",
                demux= "{demux}",
                sample_tab=sample_tab,
                config=config
        script: "../wrappers/mgi_create_samplesheet/script.py"

    rule mgi_calDemux:
        input:  sample_sheet="mgi_sample_sheet_{lane}_{demux}.txt"
        output: demultiplex_complete = "{demux}/FS2000/{lane}/demux_ready.txt",
                stats=expand("{{demux}}/FS2000/{{lane}}/{stat_filenames}.txt" \
                                            ,stat_filenames=["SequenceStat", "BarcodeStat"])

        params: tmp_dir=GLOBAL_TMPD_PATH,
                executable_file_path=GLOBAL_REF_PATH + "/general/MGI_SplitBarcode-v2/SplitBarcode-v2.0.0/linux/bin/splitBarcode",
                demux_data_dir = config["run_dir"] + "/{lane}/calFile",
                demux = "{demux}",
                lane="{lane}",
                sample_tab=sample_tab,
                run_info=expand("{run_dir}/L01/BioInfo.csv",run_dir=config["run_dir"])[0]
        threads: 60
        log: "logs/calDemux_{demux}_{lane}.log"
        # conda: "../wrappers/mgi_calDemux/env.yaml"
        script: "../wrappers/mgi_calDemux/script.py"


    # def mgi_stats_copy_input(wildcards):
    #     return expand("{demux_setting}/FS2000/{lane}/{stat_filenames}.txt" \
    #         ,lane=get_lanes_for_library(wildcards.library)
    #         ,demux_setting=get_demux_for_library(wildcards.library)
    #         ,stat_filenames = ["SequenceStat", "BarcodeStat"])
    #
    # rule mgi_stats_copy:
    #     input: mgi_stats_copy_input
    #     output: "{library}/sequencing_run_info/Stats.json",
    #     log: run="{library}/sequencing_run_info/stats_copy.log",
    #     shell:
    #         "cat {input} > {output}"

else:
    rule illumina_create_samplesheet:
        input:  run_complete_check = expand("{run_dir}/RTAComplete.txt",run_dir = config["run_dir"]),
                run_info = expand("{run_dir}/RunInfo.xml",run_dir=config["run_dir"]),
        output: samplesheet_csv = "{demux_setting}/run_samplesheet.csv",
        params: sample_tab = sample_tab,
                date = config["run_date"],
                run_name = config["run_name"],
                run_forward_read_length = config["run_forward_read_length"],
                run_reverse_read_length = config["run_reverse_read_length"],
                demux_setting= "{demux_setting}",
                run_lane_splitting = config["run_lane_splitting"]
        script: "../wrappers/create_samplesheet/script.py"

    rule illumina_bcl2fastq:
        input:  samplesheet_csv = "{demux_setting}/run_samplesheet.csv"
        output: demultiplex_complete = "{demux_setting}/demux_ready.txt",
                html  = "{demux_setting}/Stats/bcl2fastq_multiqc.html",
        params: tmp_dir = GLOBAL_TMPD_PATH,
                sample_tab = lambda wildcards: sample_tab[sample_tab['demux_setting'] == wildcards.demux_setting],
                run_lane_splitting= config["run_lane_splitting"],
                run_dir=config["run_dir"],
        log:    "logs/{demux_setting}bcl2fastq.log"
        params: library_configs = lambda wildcards: {lib_name:config["libraries"][lib_name] for lib_name in set(sample_tab[sample_tab["bcl2fastq_params_slug"] == wildcards.bcl2fastq_params_slug].library)}
        conda: "../wrappers/bcl2fastq/env.yaml"
        script: "../wrappers/bcl2fastq/script.py"


rule fastq_mv:
    input: fastq_mv_ready_input
    output: "{library}/raw_fastq/{sample_name}_R{read_num}.fastq.gz",
    params: fastq = fastq_mv_fastq_input
    threads: 60
    script: "../wrappers/fastq_mv/script.py"

rule stats_copy:
    input: fastq_mv_ready_input
    output: "{library}/sequencing_run_info/demux_info.tsv",
    params: fastq = stats_copy_input
    threads: 60
    script: "../wrappers/stats_copy/script.py"