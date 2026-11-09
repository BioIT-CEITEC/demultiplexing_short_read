######################################
# wrapper for rule: illumina_bclconvert
######################################

import os
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")

log_filename = str(snakemake.log)


def log(text):
    with open(log_filename, 'at') as f:
        f.write(text + "\n")


def run(command):
    log("## COMMAND: " + command)
    shell(command)


with open(log_filename, 'wt') as f:
    f.write("\n##\n## RULE: bclconvert \n##\n")

sample_tab = snakemake.params.sample_tab
demux_setting = sample_tab.iloc[0].to_dict()
print(demux_setting)

# BCL Convert is distributed by Illumina only as RPM/binary (not on conda)
tmp_tool_dir = os.path.join(snakemake.params.tmp_dir, "bclconvert_exec")
os.makedirs(tmp_tool_dir, exist_ok=True)
bclconvert_bin = os.path.join(tmp_tool_dir, "bcl-convert")

run("rsync -rt " + snakemake.params.bclconvert_bin + " " + bclconvert_bin + " >> " + log_filename + " 2>&1")
run("chmod +x " + bclconvert_bin)

version = str(subprocess.Popen(bclconvert_bin + " --version 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
log("## INFO: " + version)

# copy run data to tmp
run_dir = snakemake.params.run_dir
tmp_run_data = os.path.join(snakemake.params.tmp_dir, "run_tmp_data")
os.makedirs(tmp_run_data, exist_ok=True)
run("rsync -rt " + run_dir + "/* " + tmp_run_data + " >> " + log_filename + " 2>&1")

fastq_output_dir = os.path.dirname(snakemake.input.samplesheet_csv)
if fastq_output_dir == "":
    fastq_output_dir = "."

THREADS = int(snakemake.threads)
DECOMP_TH = max(1, THREADS // 3)
COMP_LVL = 9

if snakemake.params.run_lane_splitting is None:
    no_lane_splitting = " --no-lane-splitting true"
else:
    no_lane_splitting = ""

additional_options = str(demux_setting["demultiplex_additional_options"])

# --force: the output dir already exists (it holds the sample sheet)
command = bclconvert_bin \
          + " --bcl-input-directory " + tmp_run_data \
          + " --output-directory " + fastq_output_dir \
          + " --sample-sheet " + snakemake.input.samplesheet_csv \
          + " --force" \
          + no_lane_splitting \
          + " --bcl-num-conversion-threads " + str(THREADS) \
          + " --bcl-num-compression-threads " + str(THREADS) \
          + " --bcl-num-decompression-threads " + str(DECOMP_TH) \
          + " --fastq-gzip-compression-level " + str(COMP_LVL) \
          + " " + additional_options \
          + " >> " + log_filename + " 2>&1"
run(command)

with open("demux_info.tsv", 'w') as file:
    file.write(snakemake.wildcards.demux_setting + "\tall\t" + command + "\t"
               + ";".join(sample_tab['library'].unique().tolist()) + "\n")

# MultiQC's bclconvert module reads Demultiplex_Stats.csv, Quality_Metrics.csv,
# Top_Unknown_Barcodes.csv and RunInfo.xml, all written to Reports/
run("multiqc -f -z -n " + snakemake.output.html
    + " " + os.path.join(fastq_output_dir, "Reports")
    + " >> " + log_filename + " 2>&1")

run("touch " + snakemake.output.demultiplex_complete)
