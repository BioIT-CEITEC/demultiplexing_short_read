import os
from snakemake.utils import min_version
import pandas as pd

min_version("5.18.0")

GLOBAL_REF_PATH = "/mnt/ssd/ssd_3/references"

####################################
# FOLDERS
#
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

##### Config processing #####
def get_panda_sample_tab_from_config_one_lib(lib_name):
    sample_tab = pd.DataFrame.from_dict(config["libraries"][lib_name]["samples"],orient="index")
    sample_tab["library"] = lib_name
    return sample_tab


def get_panda_sample_tab_from_config(config):
    tab_list = [get_panda_sample_tab_from_config_one_lib(lib_name) for lib_name in config["libraries"].keys()]
    sample_tab = pd.concat(tab_list)
    return sample_tab

sample_tab = get_panda_sample_tab_from_config(config)

if config["run_reverse_read_length"] == 0:
    read_pair_tags = ["R1"]
else:
    read_pair_tags = ["R1","R2"]


wildcard_constraints:
    sample="[^\.]+",
    pair="R1|R2|R3|R4|SE"



####################################
# SEPARATE RULES

#include: "rules/prepare_reference.smk"
include: "rules/fastq_prepare.smk"
include: "rules/bcl2fastq.smk"

####################################
# RULE ALL
rule all:
    input:expand("{prefix}_{read_pair_tag}.fastq.gz",prefix=expand("{library}/raw_fastq/{sample}",zip,sample=sample_tab.sample_name,library=sample_tab.library),read_pair_tag=read_pair_tags),
    #input: expand(os.path.join(DIR,"raw_fastq/{sample}"))zz

