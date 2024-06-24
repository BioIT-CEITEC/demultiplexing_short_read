######################################
# wrapper for rule: bcl2fastq
######################################
import os
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")

log_filename = str(snakemake.log)

shell("conda config --add channels defaults; \
conda config --add channels bioconda; \
conda config --add channels conda-forge; \
conda config --set channel_priority strict; \
conda install -y multiqc")

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: bcl2fastq \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

sample_tab = snakemake.params.sample_tab
bcl2fastq_setting = sample_tab.iloc[0].to_dict()

print(bcl2fastq_setting)

LOAD_TH = 10
PROC_TH = 30
WRIT_TH = 10
ADAPT_STR = 0.9
MTRL = 2
MSAR = 2
COMP_LVL = 9

# base masking for special UMI processing
if len(bcl2fastq_setting["illumina_basemask"]) == 0:
  bases_mask_text = ""
else:
  bases_mask_text = " --use-bases-mask " + bcl2fastq_setting["illumina_basemask"]

if snakemake.params.run_lane_splitting == None:
  no_lane_splitting = " --no-lane-splitting"
else:
  no_lane_splitting = ""

barcode_mismatches = bcl2fastq_setting["barcode_mismatches"]
additional_options = bcl2fastq_setting["illumina_additional_options"]

run_dir = snakemake.params.run_dir

# if "config[run_sequencer_type]" == "NovaSeq":
#   bcl2fastq_args_staged_bcl_dir = os.path.join(bcl_run_dir, "Files")
# else:
#   bcl2fastq_args_staged_bcl_dir = bcl_run_dir

tmp_run_data = os.path.join(snakemake.params.tmp_dir,"run_tmp_data")
if not os.path.exists(tmp_run_data):
    os.makedirs(tmp_run_data)

command = "rsync -rt " + run_dir + "/* " + tmp_run_data + " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

fastq_output_dir = os.path.dirname(snakemake.input.samplesheet_csv)

command = "bcl2fastq -R " + tmp_run_data \
                 + " -o " + fastq_output_dir \
                 + no_lane_splitting \
                 + bases_mask_text \
                 + " --interop-dir " + tmp_run_data + "/InterOp" \
                 + " --sample-sheet " + snakemake.input.samplesheet_csv \
                 + " --loading-threads " + str(LOAD_TH) \
                 + " --processing-threads " + str(PROC_TH) \
                 + " --writing-threads " + str(WRIT_TH) \
                 + " --adapter-stringency " + str(ADAPT_STR) \
                 + " --barcode-mismatches " + str(barcode_mismatches) \
                 + " --minimum-trimmed-read-length " + str(MTRL) \
                 + " --mask-short-adapter-reads " + str(MSAR) \
                 + " " + str(additional_options) \
                 + " --fastq-compression-level " + str(COMP_LVL) \
                 + " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

if not os.path.isfile("demux_info.tsv"):
    with open("demux_info.tsv", 'w') as file:
        # Write the specified text to the file
        file.write("demux_id\tlane\trun_command\tlibraries\n")

with open("demux_info.tsv", 'w') as file:
    # Write the specified text to the file
    file.write(snakemake.wildcards.demux_setting+"\tall\t"+command+"\t"+";".join(snakemake.params.sample_tab['library'].unique().tolist())+"\n")


command = "multiqc -f -z -n "+snakemake.output.html+\
          " "+fastq_output_dir+ "/Stats/Stats.json" +\
          " >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "touch " + snakemake.output.demultiplex_complete
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# command = "multiqc -f -z -n "+snakemake.output.html+\
#           " "+snakemake.output.stats+\
#           " >> "+log_filename+" 2>&1"
# f = open(log_filename, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)
