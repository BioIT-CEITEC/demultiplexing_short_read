######################################
# wrapper for rule: bcl2fastq
######################################
import os
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: bcl2fastq \n##\n")
f.close()

library_configs = snakemake.params.library_configs

bcl2fastq_setting = list(library_configs.values())[0]

LOAD_TH = 10
PROC_TH = 30
WRIT_TH = 10
ADAPT_STR = 0.9
MTRL = 2
MSAR = 2
COMP_LVL = 4

# base masking for special UMI processing
if len(bcl2fastq_setting["base_mask_field"]) == 0:
  bases_mask_text = ""
else:
  bases_mask_text = " --use-bases-mask " + bcl2fastq_setting["base_mask_field"]

if bcl2fastq_setting["no_lane_splitting"] == True:
  no_lane_splitting = " --no-lane-splitting"
else:
  no_lane_splitting = ""

barcode_mismatches = bcl2fastq_setting["barcode_mismatches"]
additional_options = bcl2fastq_setting["additional_options"]

bcl_run_dir = os.path.dirname(snakemake.input.run_complete_check[0])

if "config[run_sequencer_type]" == "NovaSeq":
  bcl2fastq_args_staged_bcl_dir = os.path.join(bcl_run_dir, "Files")
else:
  bcl2fastq_args_staged_bcl_dir = bcl_run_dir

fastq_ouput_dir = os.path.dirname(snakemake.input.samplesheet_csv)

command = "bcl2fastq -R " + bcl2fastq_args_staged_bcl_dir \
                 + " -o " + fastq_ouput_dir \
                 + no_lane_splitting \
                 + bases_mask_text \
                 + " --interop-dir " + bcl_run_dir + "/InterOp" \
                 + " --sample-sheet " + snakemake.input.samplesheet_csv \
                 + " --loading-threads " + str(LOAD_TH) \
                 + " --processing-threads " + str(PROC_TH) \
                 + " --writing-threads " + str(WRIT_TH) \
                 + " --adapter-stringency " + str(ADAPT_STR) \
                 + " --barcode-mismatches " + str(barcode_mismatches) \
                 + " --minimum-trimmed-read-length " + str(MTRL) \
                 + " --mask-short-adapter-reads " + str(MSAR) \
                 + str(additional_options) \
                 + " --fastq-compression-level " + str(COMP_LVL) \
                 + " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
