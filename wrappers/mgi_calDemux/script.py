######################################
# wrapper for rule: bcl2fastq
######################################
import os
import pandas as pd
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")


log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: demux_MGI \n##\n")
f.close()

def get_sequencing_run_info(file_path):
    bioinfo_df = pd.read_csv(file_path)

    # Initialize variables to hold the sequence lengths
    i7_sequence_length = None
    i5_sequence_length = None

    # Loop through the DataFrame to find the relevant information
    for index, row in bioinfo_df.iterrows():
        if 'Barcode' in row.values:
            i7_sequence_length = int(row.values[1])  # Assuming the length is in the second column
        elif 'Dual Barcode' in row.values:
            i5_sequence_length = int(row.values[1])  # Assuming the length is in the second column
        elif 'Read1 Cycles' in row.values:
            read1_cycles = int(row.values[1])  # Assuming the length is in the second column
        elif 'Read2 Cycles' in row.values:
            read2_cycles = int(row.values[1])  # Assuming the length is in the second column

    return [i7_sequence_length, i5_sequence_length,read1_cycles,read2_cycles]

#run the demultiplexing
run_dir = snakemake.params.demux_data_dir
sample_tab = snakemake.params.sample_tab

tmp_run_data = os.path.join(snakemake.params.tmp_dir,"run_tmp_data",snakemake.params.lane)
if not os.path.exists(tmp_run_data):
    os.makedirs(tmp_run_data)

command = "rsync -rt " + run_dir + "/* " + tmp_run_data + " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)



# -F /var/run/user/1000/gvfs/sftp:host=147.251.158.86,user=128327/mnt/nfs/shared/000001-Instruments/MGI/C2130410230026/Raw/V350223514/L02/calFile ## path to cal folder
# -B /var/run/user/1000/gvfs/sftp:host=147.251.158.86,user=128327/mnt/share/share/710000-CEITEC/713000-cmm/713004-genomics/base/MGI/samplesheets/V350223514_L02.txt ## path to index file, concatenated i5-i7
# -r ## rev. complement both indexes
# -i 203 8 1 ## i5 index, start sequencing cycle (eg. 100+1+100+1 +1), length, errors
# -i 211 8 1 ## i5 index, start sequencing cycle (eg. i5 start + i5 length), length, errors
# -U 2 ## lane number
# -C 218 ## sequencing cycles (eg. 100+1+100+1+8+8 = 218)
# --Col 6 ## fov columns
# --Row 72 ## fov rows
# -E 3 ## cycle ends 3 = PE sequencing
# -P 111 ## length of read 1 including 1 base+
# --filter_param 2 20 20 1 1 0.75 0.75 ## default filter_mode 2 ?? 2 = on, 1 = off ?? | filter_quality_read1 20 filter_quality_read2 20 ?? min. Qscore read1, read2 ?? | filter_cyc1 110 filter_cyc2 106 ?? 1 = 100% of length ?? | filter_error_num1 28 filter_error_num2 32 ?? 0.75 = max. 25 (100-75) % errors ??; for digiMLPA --filter_param 2 20 1 1 0 0.75 0
# -o /var/run/user/1000/gvfs/sftp:host=147.251.158.86,user=128327/mnt/share/share/710000-CEITEC/713000-cmm/713004-genomics/base/MGI/V350223514/demux ## output base dir

sequencing_run_info = get_sequencing_run_info(snakemake.params.run_info)
filter_param = sample_tab[sample_tab['demux_setting'] == snakemake.params.demux].iloc[0]['MGI_filter_param']
if filter_param == "":
    filter_param = "2 20 20 1 1 0.75 0.75"
# print(filter_param)
# fastq_output_dir = os.path.dirname(snakemake.output.fastq_files[0])
# if fastq_output_dir == "":
#     fastq_output_dir = "."

command = snakemake.params.executable_file_path + " -r " \
                + " -F " + tmp_run_data \
                + " -B " + snakemake.input.sample_sheet \
                + " -i " + str(sequencing_run_info[2] + 1 + sequencing_run_info[3] + 1 + 1) \
                + " " + str(sequencing_run_info[1]) + " " + str(1) \
                + " -i " + str(sequencing_run_info[2] + 1 + sequencing_run_info[3] + 1 + 1 + sequencing_run_info[1]) \
                + " " + str(sequencing_run_info[0]) + " " + str(1) \
                + " -U " + snakemake.params.lane.replace("L0","") \
                + " -C " + str(sequencing_run_info[2] + 1 + sequencing_run_info[3] + 1 + sequencing_run_info[1] + sequencing_run_info[0]) \
                + " --Col " + "6" \
                + " --Row " + "72" \
                + " -E " + "3" \
                + " -P " + str(sequencing_run_info[3] + 1) \
                + " --filter_param " + filter_param\
                + " -o " + snakemake.params.demux\
                + " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "touch " + snakemake.output.ready_file
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
