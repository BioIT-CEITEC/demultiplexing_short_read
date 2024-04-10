from snakemake.utils import min_version
import pandas as pd
import re
import os

min_version("5.18.0")
configfile: "config.json"
GLOBAL_TMPD_PATH = config["globalTmpdPath"]
GLOBAL_REF_PATH = config["globalResources"]

#Files subdirectory check for NovaSeq outputs
if os.path.exists(config["run_dir"] + "/Files/RTAComplete.txt"):
    config["run_dir"] = config["run_dir"] + "/Files"

##### Sample table creation #####
def get_panda_sample_tab_from_config_one_lib(lib_name, run_lane_splitting_count):
    lib_config = config["libraries"][lib_name]
    sample_tab = pd.DataFrame.from_dict(lib_config["samples"],orient="index")
    sample_tab["library"] = lib_name

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
    prefix = config["run_sequencer_type"] + "_"
    for key, value in lib_config.items():
        if key.startswith(prefix):
            sample_tab[key] = value

    sample_tab["barcode_mismatches"] = lib_config["barcode_mismatches"]

    return sample_tab


def get_panda_sample_tab_from_config(config):
    run_lane_splitting_count = config.get("run_lane_splitting",None)  # Default to 4 lanes if not specified
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
    sample_tab['sample_index'] = sample_tab.groupby('demux_setting').cumcount() + 1

    # Optionally, you can drop the temporary 'combination_id' column if it's no longer needed
    sample_tab.drop(columns=['combination_id'],inplace=True)

    return sample_tab

def get_used_lanes(sample_tab,run_lane_splitting_count):

    if run_lane_splitting_count != None:
        # Extract unique libraries from sample_tab
        libraries = sample_tab['library'].unique()

        # Initialize an empty list to store the data
        libraries_lane_usage = {}

        for lib in libraries:
            # Initialize the list for lane usage for the current library
            lane_usage = []

            for lane in range(1,run_lane_splitting_count + 1):
                lane_key = f"lane{lane}"
                lane_code = f"L0{lane}"
                if sample_tab[sample_tab['library'] == lib].iloc[0][lane_key]:
                    lane_usage.append(lane_code)


            # Append the library name and its lane usage to the list
            libraries_lane_usage[lib] = lane_usage


        return libraries_lane_usage
    else:
        return None







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
    primary_lib_name = list(config["libraries"].keys())[0]
    primary_lib_raw_fastq_dir = os.path.join(primary_lib_name,"raw_fastq")
    primary_files = [f for f in os.listdir(primary_lib_raw_fastq_dir) if os.path.isfile(os.path.join(primary_lib_raw_fastq_dir,f))]

    # Function to check if all directories contain a file with a specific suffix
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

    # Remove files from primary_files based on the checks
    if not all_contain_R2:
        primary_files = [f for f in primary_files if "_R2.fastq.gz" not in f]
    if not all_contain_R3:
        primary_files = [f for f in primary_files if "_R3.fastq.gz" not in f]

    sample_tab = sample_tab.loc[sample_tab.library == primary_lib_name]
    return(sample_tab)

##### inputs to rule all #####
if "merged" in config and config["merged"]:
    sample_tab = get_sample_tab_for_merge(config)
    # TODO:
    all_sample_fastq_files = [os.path.join(library_names,\
                                      "qc_reports",\
                                      re.sub("_R..fastq.gz$","",f),\
                                      "/raw_fastqc/",\
                                      re.sub(".fastq.gz$","_fastqc.html",f)) for f in primary_files]


    library_output = list(config["library_output"].keys())[0]
    library_names = library_output

else:
    sample_tab = get_panda_sample_tab_from_config(config)

    # sample_tab = sample_tab.iloc[:15]
    # print(sample_tab)

    sample_file_tab = sample_tab.reindex(sample_tab.index.repeat(sample_tab['read_output_count'])) \
        .assign(read_num=lambda x: x.groupby(['library', 'sample_name']).cumcount() + 1) \
        .reset_index(drop=True)

    per_library_used_lanes = get_used_lanes(sample_file_tab,config.get("run_lane_splitting",None))
    library_names = set(sample_tab["library"])
    resulting_fastq_files = expand("{library}/raw_fastq/{sample_name}_R{read_num}.fastq.gz",zip \
        ,library=sample_file_tab.library \
        ,sample_name=sample_file_tab.sample_name \
        ,read_num=sample_file_tab.read_num)

    library_output = "-"
    # print(per_library_used_lanes)

##### wildcard_constraints #####
wildcard_constraints:
    sample = "|".join(set(sample_tab.sample_name.tolist())),
    library = "|".join(set(sample_tab.library.tolist())),
    bcl2fastq_params_slug = "bcl2fastqslug_[a-zA-Z0-9_-]*",
    demux="demux_[0-9]"


rule all:
    input: fastq_files = resulting_fastq_files,
           stats = expand("{library}/sequencing_run_info/Stats.json",library = library_names)

##### Modules #####

include: "rules/merge.smk"
include: "rules/demultiplexing.smk"
include: "rules/fastqc.smk"
# include: "rules/check_adaptors.smk"

