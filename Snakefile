from snakemake.utils import min_version
import pandas as pd
import re
import os

min_version("5.18.0")
configfile: "config.json"
GLOBAL_TMPD_PATH = config["globalTmpdPath"]

#Files subdirectory check for NovaSeq outputs
if os.path.exists(config["run_dir"] + "/Files/RTAComplete.txt"):
    config["run_dir"] = config["run_dir"] + "/Files"

##### Sample table creation #####
def get_panda_sample_tab_from_config_one_lib(lib_name):
    sample_tab = pd.DataFrame.from_dict(config["libraries"][lib_name]["samples"],orient="index")
    sample_tab["library"] = lib_name
    bcl2fastq_params_slug = ("bcl2fastqslug_" + str(config["libraries"][lib_name]["barcode_mismatches"]) + "" + \
                            str(config["libraries"][lib_name]["additional_options"]) + "" + \
                            str(config["libraries"][lib_name]["R1FastQMask"]) + "" + \
                            str(config["libraries"][lib_name]["R2FastQMask"]) + "" + \
                            str(config["libraries"][lib_name]["I1Mask"]) + "" + \
                            str(config["libraries"][lib_name]["I2Mask"]) + "" + \
                            str(config["libraries"][lib_name]["UmiMask"]))
    sample_tab["bcl2fastq_params_slug"] = re.sub("[^a-zA-Z0-9_-]","_",bcl2fastq_params_slug)
    return sample_tab

def get_panda_sample_tab_from_config(config):
    tab_list = [get_panda_sample_tab_from_config_one_lib(lib_name) for lib_name in config["libraries"].keys()]
    sample_tab = pd.concat(tab_list)
    return sample_tab

sample_tab = get_panda_sample_tab_from_config(config)
sample_tab['read_output_count'] = sample_tab.apply(lambda row: 1 if config["run_reverse_read_length"] == 0 else 2, axis=1)
for lib_name in config["libraries"].keys():
    # Check if base_mask_field is not empty
    if config["libraries"][lib_name]["R2FastQMask"] != "":
        # Count the 'y' characters in base_mask_field
        y_count = config["libraries"][lib_name]["R2FastQMask"].upper().count('Y')
        if y_count != 0:
            # Update the DataFrame
            sample_tab.loc[sample_tab['library'] == lib_name, 'read_output_count'] = 2
        else:
            sample_tab.loc[sample_tab['library'] == lib_name, 'read_output_count'] = 1
    y_count = config["libraries"][lib_name]["UmiMask"].upper().count('Y')
    if y_count != 0:
        # Update the DataFrame
        sample_tab.loc[sample_tab['library'] == lib_name, 'read_output_count'] += 1

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

    all_merge_sample_fastqc_files = [os.path.join(library_output,\
                                      "qc_reports",\
                                      re.sub("_R..fastq.gz$","",f),\
                                      "/raw_fastqc/",\
                                      re.sub(".fastq.gz$","_fastqc.html",f)) for f in primary_files]
    library_names = library_output
    sample_tab = sample_tab.loc[sample_tab.library == primary_lib_name]
else:
    sample_file_tab = sample_tab.reindex(sample_tab.index.repeat(sample_tab['read_output_count'])) \
        .assign(read_num=lambda x: x.groupby(['library', 'sample_name']).cumcount() + 1) \
        .reset_index(drop=True)

    library_names = set(config["libraries"].keys())

rule all:
    input: fastq_files = expand("{library}/qc_reports/raw_fastq_multiqc.html",library = library_names),
           stats = expand("{library}/sequencing_run_info/Stats.json",library = library_names)

##### Modules #####

include: "rules/merge.smk"
include: "rules/demultiplexing.smk"
include: "rules/fastqc.smk"
# include: "rules/check_adaptors.smk"

