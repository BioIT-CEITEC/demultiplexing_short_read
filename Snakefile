from snakemake.utils import min_version
import pandas as pd
import re
import os

min_version("5.18.0")
configfile: "config.json"

#Files subdirectory check for NovaSeq outputs
if os.path.exists(config["run_dir"] + "/Files/RTAComplete.txt"):
    config["run_dir"] = config["run_dir"] + "/Files"

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
sample_tab = sample_tab.set_index(pd.RangeIndex(start=1,stop=len(sample_tab.index) + 1))

#duplicated sample name check
existing_name_check_dict = {}

for index, row in sample_tab.iterrows():
    if row["sample_name"] in existing_name_check_dict.keys():
        sample_tab.loc[sample_tab['library'] == row["library"], 'bcl2fastq_params_slug'] = row["bcl2fastq_params_slug"] + "_" + str(index)
    existing_name_check_dict[row["sample_name"]] = True

sample_tab["slug_id"] = sample_tab.groupby("bcl2fastq_params_slug").transform(lambda x: range(1,len(x.index) + 1))["sample_name"]

if "library_output" in config:
    library_output = list(config["library_output"].keys())[0]
else:
    library_output = "-"

##### wildcard_constraints #####
wildcard_constraints:
    sample = "|".join(set(sample_tab.sample_name.tolist())),
    library = "|".join(set(sample_tab.library.tolist())),
    bcl2fastq_params_slug = "bcl2fastqslug_[a-zA-Z0-9_-]*"

##### inputs to rule all #####
if "merged" in config and config["merged"]:
    primary_lib_name = list(config["libraries"].keys())[0]
    primary_lib_raw_fastq_dir = os.path.join(primary_lib_name,"raw_fastq")
    primary_files = [f for f in os.listdir(primary_lib_raw_fastq_dir) if os.path.isfile(os.path.join(primary_lib_raw_fastq_dir,f))]
    all_sample_inputs = [os.path.join(library_output,"raw_fastq",f) for f in primary_files]
    library_names = library_output
    sample_tab = sample_tab.loc[sample_tab.library == primary_lib_name]
else:

    ##### All resulting fastqs #####
    single_end_samples = [row["library"] + "/raw_fastq/" + row["sample_name"] + ".fastq.gz" for  index, row in sample_tab.iterrows() if config["libraries"][row["library"]]["lib_reverse_read_length"] == 0]
    first_pair_end_samples = [row["library"] + "/raw_fastq/" + row["sample_name"] + "_R1.fastq.gz" for  index, row in sample_tab.iterrows() if config["libraries"][row["library"]]["lib_reverse_read_length"] != 0]
    second_pair_end_samples = [row["library"] + "/raw_fastq/" + row["sample_name"] + "_R2.fastq.gz" for  index, row in sample_tab.iterrows() if config["libraries"][row["library"]]["lib_reverse_read_length"] != 0]
    all_sample_inputs = single_end_samples + first_pair_end_samples + second_pair_end_samples
    library_names = set(config["libraries"].keys())

rule all:
    input: fastq_files = all_sample_inputs,
           stats = expand("{library_name}/sequencing_run_info/Stats.json",library_name = library_names)

##### Modules #####

include: "rules/merge.smk"
include: "rules/fastq_prepare.smk"
include: "rules/bcl2fastq.smk"





