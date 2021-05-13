import os
from snakemake.utils import min_version
import pandas as pd

min_version("5.18.0")

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
#DIR = os.path.join("/mnt/ssd/ssd_1/snakemake" ,config["libraries"]["vojta2"]["samples"])

#sample_tab = os.path.join()

if config["run_reverse_read_length"] == 0:
    read_pair_tags = ["R1"]
else:
    read_pair_tags = ["R1","R2"]

##### Target rules #####

rule all:
    input:expand("{prefix}_{read_pair_tag}.fastq.gz",prefix=expand("{library}/raw_fastq/{sample}",zip,sample=sample_tab.sample_name,library=sample_tab.library),read_pair_tag=read_pair_tags),
    #input: expand(os.path.join(DIR,"raw_fastq/{sample}"))zz


##### Modules #####

include: "rules/bcl2fastq.smk"


