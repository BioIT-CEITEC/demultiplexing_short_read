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

bcl2fastq_params_slug_list = [str(one_library_config["barcode_mismatches"]) + "__" + \
                                       str(one_library_config["additional_options"]) + "__" + \
                                       str(one_library_config["base_mask_field"]) + "__" + \
                                       str(one_library_config["no_lane_splitting"]) for one_library_config in library_configs.values()]

if len(set(bcl2fastq_params_slug_list)) == 1:
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

    fastq_ouput_dir = os.path.dirname(snakemake.output.fastqs[0])

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

    for lib_name in library_configs.keys():
        command = "mkdir -p " + lib_name
        f = open(log_filename, 'at')
        f.write("## COMMAND: " + command + "\n")
        f.close()
        shell(command)

        command = "cp -r " + fastq_ouput_dir + "/Stats " + lib_name + "/sequencing_run_info/"
        f = open(log_filename, 'at')
        f.write("## COMMAND: " + command + "\n")
        f.close()
        shell(command)

        command = "cp -r " + fastq_ouput_dir + "/Reports " + lib_name + "/sequencing_run_info/"
        f = open(log_filename, 'at')
        f.write("## COMMAND: " + command + "\n")
        f.close()
        shell(command)

    #if snakemake.params.run_log_dir:
    #    read_sum = 0
    #    for fastqc_file in snakemake.input.html:
    #        f = open(fastqc_file, 'r')
    #        for line in f:
    #            if "Total Sequences</td><td>" in line:
    #                num = re.sub('.*Total Sequences</td><td>([0-9]*)</td>.*',r'\1',line)
    #                num = int(num)
    #        read_sum = read_sum + num

        #f = open(snakemake.params.run_log_dir + "/per_library_read_sum.tsv", 'a+')
        #f.write("\t"+str(read_sum)+"\n")
        #f.close()

#
# if [ $? -eq 0 ]
# then
#   cp -rt ${BCL2FASTQ_LOGS_OUTPUT}/${RUNDIRNAME} ${STAGED_DIR}/raw_fastq_temp/${RUNDIRNAME}/*/
#   LOG_DIR=$(dirname "${LOG}")
#   echo "## COMMAND: cp -r ${LOG_DIR} ${BCL2FASTQ_LOGS_OUTPUT}"
#   echo "## COMMAND: cp -rt ${BCL2FASTQ_LOGS_OUTPUT}/${RUNDIRNAME} ${STAGED_DIR}/raw_fastq_temp/${RUNDIRNAME}/*/"
#   cp -r ${LOG_DIR} ${BCL2FASTQ_LOGS_OUTPUT}
#   cp -rt ${BCL2FASTQ_LOGS_OUTPUT}/${RUNDIRNAME} ${STAGED_DIR}/raw_fastq_temp/${RUNDIRNAME}/*/
#   ORIG_RUN_NAME=$(basename "${ORIGINAL_BCL_DIR}")
#   rm -rf ${TEMP_RAW_BCL_DIR}${ORIG_RUN_NAME}
#   rm -rf ${LOG_DIR}
#   touch ${STAGED_DIR}/raw_fastq_temp/${RUNDIRNAME}/task_finished_OK
# else
#   rm -f $STLOG
#   cp $LOG $STLOG
# fi
#
# echo "## COMMAND: conda deactivate"
# conda deactivate


# conda deactivate
