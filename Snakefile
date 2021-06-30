from snakemake.utils import min_version
import pandas as pd

min_version("5.18.0")

##### Sample table creation #####
def get_panda_sample_tab_from_config_one_lib(lib_name):
    sample_tab = pd.DataFrame.from_dict(config["libraries"][lib_name]["samples"],orient="index")
    sample_tab["library"] = lib_name
    bcl2fastq_params_slug = "bcl2fastqslug_" + str(config["libraries"][lib_name]["barcode_mismatches"]) + "" + \
                                  str(config["libraries"][lib_name]["additional_options"]) + "" + \
                                  str(config["libraries"][lib_name]["base_mask_field"]) + "" + \
                                  str(config["libraries"][lib_name]["no_lane_splitting"])
    sample_tab["bcl2fastq_params_slug"] = re.sub("[^a-zA-Z0-9_-]","_",bcl2fastq_params_slug)
    return sample_tab

def get_panda_sample_tab_from_config(config):
    tab_list = [get_panda_sample_tab_from_config_one_lib(lib_name) for lib_name in config["libraries"].keys()]
    sample_tab = pd.concat(tab_list)
    return sample_tab

sample_tab = get_panda_sample_tab_from_config(config)
sample_tab = sample_tab.set_index(pd.RangeIndex(start = 1,stop = len(sample_tab.index)+1))
sample_tab["slug_id"] = sample_tab.groupby("bcl2fastq_params_slug").transform(lambda x: range(1,len(x.index)+1))["sample_name"]

##### All resulting fastqs #####
single_end_samples = [row["library"] + "/raw_fastq/" + row["sample_name"] + ".fastq.gz" for  index, row in sample_tab.iterrows() if config["libraries"][row["library"]]["lib_reverse_read_length"] == 0]
first_pair_end_samples = [row["library"] + "/raw_fastq/" + row["sample_name"] + "_R1.fastq.gz" for  index, row in sample_tab.iterrows() if config["libraries"][row["library"]]["lib_reverse_read_length"] != 0]
second_pair_end_samples = [row["library"] + "/raw_fastq/" + row["sample_name"] + "_R2.fastq.gz" for  index, row in sample_tab.iterrows() if config["libraries"][row["library"]]["lib_reverse_read_length"] != 0]
all_sample_inputs = single_end_samples + first_pair_end_samples + second_pair_end_samples

##### wildcard_constraints #####
wildcard_constraints:
    sample = "|".join(set(sample_tab.sample_name.tolist())),
    library = "|".join(set(sample_tab.library.tolist())),
    bcl2fastq_params_slug = "bcl2fastqslug_[a-zA-Z0-9_-]*",

rule all:
    input: all_sample_inputs

##### Rules #####

rule create_samplesheet:
    params: sample_tab = lambda wildcards: sample_tab[sample_tab["bcl2fastq_params_slug"] == wildcards.bcl2fastq_params_slug],
            date = config["run_date"],
            run_name = config["run_name"],
            run_forward_read_length = config["run_forward_read_length"],
            run_reverse_read_length = config["run_reverse_read_length"]
    output: samplesheet_csv = config["run_name"] + "/{bcl2fastq_params_slug}/run_samplesheet.csv",
    script: "wrappers/create_samplesheet/script.py"


rule bcl2fastq:
    input:  run_complete_check = expand("{run_dir}/RTAComplete.txt",run_dir = config["run_dir"]),
            samplesheet_csv= config["run_name"] + "/{bcl2fastq_params_slug}/run_samplesheet.csv"
    output: demultiplex_complete_check = config["run_name"] + "/{bcl2fastq_params_slug}/Reports/html/index.html",
    params: library_configs = lambda wildcards: {lib_name:config["libraries"][lib_name] for lib_name in set(sample_tab[sample_tab["bcl2fastq_params_slug"] == wildcards.bcl2fastq_params_slug].library)}
    log:    config["run_name"] + "/{bcl2fastq_params_slug}/bcl2fastq.log"
    conda: "wrappers/bcl2fastq/env.yaml"
    script: "wrappers/bcl2fastq/script.py"


rule fastq_mv:
    input:  demultiplex_complete_check = expand(config["run_name"] + "/{bcl2fastq_params_slug}/Reports/html/index.html",bcl2fastq_params_slug = list(pd.unique(sample_tab['bcl2fastq_params_slug']))),
    output: fastqs_out = expand(config["run_name"] + "/{sample}_S{sample_index}_R1_001.fastq.gz",zip,sample = sample_tab.sample_name\
                                                                            ,sample_index = range(1,len(sample_tab.index)+1))
    params: fastqs_in = expand(config["run_name"] + "/{bcl2fastq_params_slug}/{sample}_S{sample_index}_R1_001.fastq.gz",zip\
                                                                            ,sample = sample_tab.sample_name\
                                                                            ,sample_index = sample_tab.slug_id\
                                                                            ,bcl2fastq_params_slug = sample_tab.bcl2fastq_params_slug)
    script: "wrappers/fastq_mv/script.py"

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


