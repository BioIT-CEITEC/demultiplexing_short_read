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

def get_sample_ID_for_library(library,sample_name):
    sample_ID = sample_tab.loc[(sample_tab['library'] == library) & (sample_tab['sample_name'] == sample_name)].iloc[0]['sample_ID']
    return sample_ID

def fastq_mv_ready_input(wildcards):
    if config["run_sequencer_type"] == "AVITI":
        input_list = expand("{demux_setting}/demux_ready.txt" \
            ,demux_setting=get_demux_for_library(wildcards.library))
    elif config["run_sequencer_type"] == "MGI":
        input_list = expand("{demux_setting}/{run_dir_id}/{lane}/demux_ready.txt" \
            ,lane=get_lanes_for_library(wildcards.library)
            ,demux_setting=get_demux_for_library(wildcards.library)
            ,run_dir_id = config["run_dir"])
    else:
        input_list = expand("{demux_setting}/demux_ready.txt" \
            ,demux_setting=get_demux_for_library(wildcards.library))

    return input_list


def fastq_mv_fastq_input(wildcards):
    if config["run_sequencer_type"] == "AVITI":
        input_list = expand("{demux_setting}/Samples/"+wildcards.library+"/"+wildcards.sample_name+"___{sample_ID}/"+wildcards.sample_name+"___{sample_ID}_R"+wildcards.read_num+".fastq.gz"\
            ,demux_setting=get_demux_for_library(wildcards.library)
            ,sample_ID=get_sample_ID_for_library(wildcards.library,wildcards.sample_name))
    elif config["run_sequencer_type"] == "MGI":
        input_list = expand("{demux_setting}/{run_dir_id}/{lane}/{run_dir_id}_{lane}_" + wildcards.sample_name + "___{sample_ID}_" + wildcards.read_num + ".fq.gz" \
            ,lane=get_lanes_for_library(wildcards.library)
            ,demux_setting=get_demux_for_library(wildcards.library)
            ,sample_ID=get_sample_ID_for_library(wildcards.library,wildcards.sample_name)
            ,run_dir_id = config["run_dir"])
    else:
        input_list = expand("{demux_setting}/" + wildcards.sample_name + "___{sample_ID}_S{sample_index}{lane}_R" + wildcards.read_num + "_001.fastq.gz" \
            ,lane=get_lanes_for_library(wildcards.library)
            ,demux_setting=get_demux_for_library(wildcards.library)
            ,sample_ID=get_sample_ID_for_library(wildcards.library,wildcards.sample_name)
            ,sample_index=get_sample_index_for_library(wildcards.library,wildcards.sample_name))
    return input_list

def stats_copy_input(wildcards):
    if config["run_sequencer_type"] == "AVITI":
        input_list = expand("{demux_setting}/{stat_filenames}" \
            ,demux_setting=get_demux_for_library(wildcards.library)
            ,stat_filenames=["*.html", "*.json","*.csv","*.json"])
    elif config["run_sequencer_type"] == "MGI":
        input_list = expand(["{demux_setting}/{run_dir_id}/{lane}/SequenceStat.txt",
         "{demux_setting}/{run_dir_id}/{lane}/BarcodeStat.txt",
         "mgi_sample_sheet_{lane}_{demux_setting}.txt"],
            zip,
            demux_setting=get_demux_for_library(wildcards.library),
            lane=get_lanes_for_library(wildcards.library),
            run_dir_id = config["run_dir"])
    else:
        input_list = expand("{demux_setting}/{stat_filenames}" \
            ,demux_setting=get_demux_for_library(wildcards.library)
            ,stat_filenames=["Reports", "Stats","run_samplesheet.csv","Stats"])

    return input_list

def nread_file_input(wildcards):
    if config["run_sequencer_type"] == "AVITI":
        input_list = expand("{demux_setting}/IndexAssignment.csv" \
            ,demux_setting=get_demux_for_library(wildcards.library))
    elif config["run_sequencer_type"] == "MGI":
        input_list = expand("{demux_setting}/{run_dir_id}/{lane}/BarcodeStat.txt",zip \
            ,lane=get_lanes_for_library(wildcards.library)
            ,demux_setting=get_demux_for_library(wildcards.library)
            ,run_dir_id = config["run_dir"])
    else:
        input_list = expand("{demux_setting}/Stats/DemultiplexingStats.xml" \
            ,demux_setting=get_demux_for_library(wildcards.library))
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
        threads: 20
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
        output: demultiplex_complete = expand("{{demux}}/{run_dir_id}/{{lane}}/demux_ready.txt" \
                                            ,run_dir_id = config["run_dir"])[0],
                stats = expand("{{demux}}/{run_dir_id}/{{lane}}/{stat_filenames}.txt" \
                                            ,stat_filenames=["SequenceStat", "BarcodeStat"],run_dir_id = config["run_dir"])

        params: tmp_dir=GLOBAL_TMPD_PATH,
                executable_file_path=GLOBAL_REF_PATH + "/general/MGI_SplitBarcode-v2/SplitBarcode-v2.0.0/linux/bin/splitBarcode",
                run_dir_id = config["run_dir"],
                demux_data_dir = config["run_dir"] + "/{lane}/calFile",
                demux = "{demux}",
                lane="{lane}",
                sample_tab=sample_tab,
                run_info=expand("{run_dir}/L01/BioInfo.csv",run_dir=config["run_dir"])[0],
                config=config
        threads: 60
        log: "logs/calDemux_{demux}_{lane}.log"
        # conda: "../wrappers/mgi_calDemux/env.yaml"
        script: "../wrappers/mgi_calDemux/script.py"

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
    output: fastq = "{library}/raw_fastq/{sample_name}_R{read_num}.fastq.gz",
    params: fastq = fastq_mv_fastq_input,
            run_sequencer_type = config["run_sequencer_type"],
            date_id = DEMUX_DATE_ID
    threads: 1
    script: "../wrappers/fastq_mv/script.py"

rule stats_copy:
    input: in_stats_file = fastq_mv_ready_input
    output: nread_json = "{library}/sequencing_run_info/samplesNumberReads.json",
    params: stats_files = stats_copy_input,
            nread_file = nread_file_input,
            sample_tab = lambda wildcards: sample_tab[sample_tab["library"] == wildcards.library],
            sequencer_type = config["run_sequencer_type"]
    threads: 60
    script: "../wrappers/stats_copy/script.py"
