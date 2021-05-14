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
print(sample_tab)
#sample_tab = os.path.join()

if config["run_reverse_read_length"] == 0:
    read_pair_tags = [""]
else:
    read_pair_tags = ["_R1","_R2"]

##### Target rules #####

rule all:
    input:expand("{prefix}{read_pair_tag}.fastq.gz",prefix=expand("{library}/raw_fastq/{sample}",zip,sample=sample_tab.sample_name,library=sample_tab.library),read_pair_tag=read_pair_tags),

##### Rules #####

rule create_samplesheet:
    params: sample_tab = sample_tab,
            date = config["run_date"],
            run_name = config["run_name"],
            run_forward_read_length = config["run_forward_read_length"],
            run_reverse_read_length = config["run_reverse_read_length"]
    output: samplesheet_csv = "run_samplesheet.csv",
    script: "wrappers/create_samplesheet/script.py"

rule bcl2fastq:
    input:  run_complete_check = expand("{run_dir}/RTAComplete.txt",run_dir = config["run_dir"]),
            samplesheet_csv= "run_samplesheet.csv"
    output: fastqs = expand(config["run_name"] + "/{sample}_S{sample_index}_R1_001.fastq.gz",zip,sample = sample_tab.sample_name\
                                                                            ,sample_index = range(1,len(sample_tab.index)+1))
    params: library_configs = config["libraries"]
    log:    "bcl2fastq.log"
    conda: "wrappers/bcl2fastq/env.yaml"
    script: "wrappers/bcl2fastq/script.py"

rule fastq_prepare_SE:
    input:  in_filename = lambda wildcards: expand("{run_name}/{sample}_S{sample_index}_R1_001.fastq.gz",run_name = config["run_name"]\
                    ,sample = wildcards.sample\
                    ,sample_index = sample_tab.loc[(sample_tab["sample_name"] == wildcards.sample) & (sample_tab.library == wildcards.library)].index[0])
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
                    ,sample_index = sample_tab.loc[(sample_tab["sample_name"] == wildcards.sample) & (sample_tab.library == wildcards.library)].index[0])
    output: R1 = "{library}/raw_fastq/{sample}_R1.fastq.gz",
            R2 = "{library}/raw_fastq/{sample}_R2.fastq.gz"
    log:    "{library}/sample_logs/{sample}/fastq_prepare_PE.log"
    params: umi = lambda wildcards: config["libraries"][wildcards.library]["UMI"],
            run_name = config["run_name"],
    threads:  1
    conda:  "wrappers/fastq_prepare_PE/env.yaml"
    script: "wrappers/fastq_prepare_PE/script.py"


