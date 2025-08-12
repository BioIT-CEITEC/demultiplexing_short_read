from snakemake.utils import min_version
import pandas as pd
import re
import os
from datetime import datetime

min_version("5.18.0")
configfile: "config.json"
GLOBAL_TMPD_PATH = config["globalTmpdPath"]
GLOBAL_REF_PATH = config["globalResources"]
DEMUX_DATE_ID = datetime.now().strftime('%Y%m%d%H%M')

#Files subdirectory check for NovaSeq outputs
if os.path.exists(config["run_dir"] + "/Files/RTAComplete.txt"):
    config["run_dir"] = config["run_dir"] + "/Files"

##### Sample table creation #####
def get_panda_sample_tab_from_config_one_lib(lib_name, run_lane_splitting_count):
    lib_config = config["libraries"][lib_name]
    sample_tab = pd.DataFrame.from_dict(lib_config["samples"],orient="index")
    sample_tab["library"] = lib_name

    sample_tab['sample_ID'] = sample_tab.index.astype(str)
    sample_tab['sample_name_full'] = sample_tab['sample_name'] + '___' + sample_tab['sample_ID']

    # add library size
    sample_tab['lib_forward_read_length'] = lib_config["lib_forward_read_length"]
    sample_tab['lib_reverse_read_length'] = lib_config["lib_reverse_read_length"]
    sample_tab['cut_reads_to_lib_size'] = lib_config["cut_reads_to_lib_size"]

    # Adding the read output counts column
    sample_tab['read_output_count'] = sample_tab.apply(lambda row: 1 if config[
                                                                            "run_reverse_read_length"] == 0 else 2,axis=1)
    if config["run_sequencer_type"] == "AVITI":
        # Check if base_mask_field is not empty
        if config["libraries"][lib_name]["AVITI_R2FastQMask"] != "":
            # Count the 'y' characters in base_mask_field
            y_count = config["libraries"][lib_name]["AVITI_R2FastQMask"].upper().count('Y')
            if y_count != 0:
                # Update the DataFrame
                sample_tab.loc[sample_tab['library'] == lib_name, 'read_output_count'] = 2
            else:
                sample_tab.loc[sample_tab['library'] == lib_name, 'read_output_count'] = 1
        y_count = config["libraries"][lib_name]["AVITI_UmiMask"].upper().count('Y')
        if y_count != 0:
            # Update the DataFrame
            sample_tab.loc[sample_tab['library'] == lib_name, 'read_output_count'] += 1

    # Adding the lane columns
    if run_lane_splitting_count != None:
        for lane in range(1,run_lane_splitting_count + 1):
            lane_key = f"lane{lane}"
            sample_tab[lane_key] = lib_config["lib_lane_splitting"].get(lane_key,True)

    # Adding the demultiplexing settings columns
    if config["run_sequencer_type"] == "AVITI" or config["run_sequencer_type"] == "MGI":
        prefix = config["run_sequencer_type"] + "_"
    else:
        prefix = "illumina_"

    for key, value in lib_config.items():
        if key.startswith(prefix):
            sample_tab[key] = value

    sample_tab["barcode_mismatches"] = lib_config["barcode_mismatches"]
    sample_tab["demultiplex_additional_options"] = lib_config["demultiplex_additional_options"]

    return sample_tab


def get_panda_sample_tab_from_config(config):
    run_lane_splitting_count = config.get("run_lane_splitting",None)
    tab_list = [get_panda_sample_tab_from_config_one_lib(lib_name,run_lane_splitting_count) for lib_name in config["libraries"].keys()]
    sample_tab = pd.concat(tab_list,ignore_index=True)


    # Identify columns that start with the run_sequencer_type prefix
    prefix_columns = [col for col in sample_tab.columns if col.startswith(config["run_sequencer_type"] + "_")]
    prefix_columns.append("barcode_mismatches")

    # Create a unique identifier for each combination of values in the prefix_columns
    sample_tab['combination_id'] = sample_tab[prefix_columns].astype(str).agg('-'.join,axis=1)

    # Map each unique combination to a unique demux_setting value
    combination_to_demux = {combination: f'demux_{i + 1}' for i, combination in
                            enumerate(sample_tab['combination_id'].unique())}
    sample_tab['demux_setting'] = sample_tab['combination_id'].map(combination_to_demux)

    # Sort by lanes (lane1, lane2, ..., laneN) in descending order so that rows with lane1=True appear first, then lane2=True, etc.
    if run_lane_splitting_count != None:
        lane_columns = [f"lane{lane}" for lane in range(1, run_lane_splitting_count + 1)]
        sample_tab.sort_values(by=lane_columns, ascending=False, inplace=True)

    sample_tab['sample_index'] = sample_tab.groupby('demux_setting').cumcount() + 1

    # Optionally, you can drop the temporary 'combination_id' column if it's no longer needed
    sample_tab.drop(columns=['combination_id'],inplace=True)

    return sample_tab

# sample_tab = sample_tab.set_index(pd.RangeIndex(start=1,stop=len(sample_tab.index) + 1))
#
# #duplicated sample name check
# existing_name_check_dict = {}
#
# for index, row in sample_tab.iterrows():
#     if row["sample_name"] in existing_name_check_dict.keys():
#         sample_tab.loc[sample_tab['library'] == row["library"], 'bcl2fastq_params_slug'] = row["bcl2fastq_params_slug"] + "_" + str(index)
#     existing_name_check_dict[row["sample_name"]] = True
#
# sample_tab["slug_id"] = sample_tab.groupby("bcl2fastq_params_slug").transform(lambda x: range(1,len(x.index) + 1))["sample_name"]

def get_sample_tab_for_merge(config):
    merged_lib_name = list(config["library_output"].keys())[0]
    lib_config = config["library_output"][merged_lib_name]
    sample_tab = pd.DataFrame.from_dict(lib_config["samples"],orient="index")
    sample_tab["library"] = merged_lib_name
    sample_tab['sample_ID'] = sample_tab.index.astype(str)

    def all_dirs_contain_suffix(suffix):
        for lib_name in config["libraries"].keys():
            lib_raw_fastq_dir = os.path.join(lib_name,"raw_fastq")
            if not any(suffix in f for f in os.listdir(lib_raw_fastq_dir) if
                       os.path.isfile(os.path.join(lib_raw_fastq_dir,f))):
                return False
        return True

    # Check for each file type independently
    all_contain_R2 = all_dirs_contain_suffix("_R2.fastq.gz")
    all_contain_R3 = all_dirs_contain_suffix("_R3.fastq.gz")

    sample_tab['read_output_count'] = 1
    if all_contain_R2:
        sample_tab['read_output_count'] = 2
    if all_contain_R3:
        sample_tab['read_output_count'] = 3

    return sample_tab


##### inputs to rule all #####
if "merged" in config and config["merged"]:
    sample_tab = get_sample_tab_for_merge(config)

    sample_file_tab = sample_tab.reindex(sample_tab.index.repeat(sample_tab['read_output_count'])) \
        .assign(read_num=lambda x: x.groupby(['library', 'sample_name']).cumcount() + 1) \
        .reset_index(drop=True)

    library_output = list(config["library_output"].keys())[0]

    resulting_fastq_files = expand("{library}/raw_fastq/{sample_name}_R{read_num}.fastq.gz",zip \
        ,library=sample_file_tab.library \
        ,sample_name=sample_file_tab.sample_name \
        ,read_num=sample_file_tab.read_num)

    library_names = library_output

else:
    sample_tab = get_panda_sample_tab_from_config(config)
    sample_tab.to_csv('sample_tab.csv', index=False)

    sample_file_tab = sample_tab.reindex(sample_tab.index.repeat(sample_tab['read_output_count'])) \
        .assign(read_num=lambda x: x.groupby(['library', 'sample_name']).cumcount() + 1) \
        .reset_index(drop=True)

    library_names = set(sample_tab["library"])
    resulting_fastq_files = expand("{library}/raw_fastq/{sample_name}_R{read_num}.fastq.gz",zip \
        ,library=sample_file_tab.library \
        ,sample_name=sample_file_tab.sample_name \
        ,read_num=sample_file_tab.read_num)

    library_output = "-"

##### wildcard_constraints #####
wildcard_constraints:
    sample = "|".join(set(sample_tab.sample_name.tolist())),
    library = "|".join(set(sample_tab.library.tolist())),
    bcl2fastq_params_slug = "bcl2fastqslug_[a-zA-Z0-9_-]*",
    demux="demux_[0-9]"


rule all:
    input: fastq_files = resulting_fastq_files,
           stats = expand("{library}/sequencing_run_info/samplesNumberReads.json",library = library_names)

##### Modules #####

if "merged" in config and config["merged"]:
    include: "rules/merge.smk"
else:
    include: "rules/demultiplexing.smk"
include: "rules/fastqc.smk"
# include: "rules/check_adaptors.smk"

