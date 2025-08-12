#############################################################
# wrapper for rule: fastq_mv
#############################################################
import os
import re
from snakemake.shell import shell
from datetime import datetime

shell.executable("/bin/bash")

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: fastq_mv \n##\n")
f.close()

#create output dir

command = "mkdir -p " + os.path.dirname(snakemake.output.fastq)
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)



in_fastq_list = []
for in_file in snakemake.params.fastq:
    if os.path.isfile(in_file):
        in_fastq_list.append(in_file)

if len(in_fastq_list) == 0:
    command = "touch " + snakemake.output.fastq
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
else:
    # Determine target read length based on read number
    # Extract read number from output filename
    read_num_match = re.search(r'_R(\d+)\.fastq\.gz$', snakemake.output.fastq)
    if read_num_match:
        read_num = int(read_num_match.group(1))
    else:
        raise ValueError(f"Could not extract read number from filename: {snakemake.output.fastq}")
    
    # Get the appropriate read length
    if read_num == 1:
        target_length = snakemake.params.lib_forward_read_length
    elif read_num == 2:
        target_length = snakemake.params.lib_reverse_read_length
    else:
        # For read 3 or higher (e.g., UMI reads), do not set any size, so it won't be trimmed
        target_length = None
    
    # Concatenate input files if multiple
    if len(in_fastq_list) == 1:
        concatenated_fastq = in_fastq_list[0]
    else:
        concatenated_fastq = snakemake.output.fastq + ".tmp_concat"
        command = "cat " + " ".join(in_fastq_list) + " > " + concatenated_fastq
        f = open(log_filename, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)
    
    # If target_length is 0 or None, just copy/move the file without trimming
    if target_length == 0 or target_length is None:
        command = "mv " + concatenated_fastq + " " + snakemake.output.fastq
        f = open(log_filename, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)
    else:
        # Use cutadapt to trim reads to target length
        cutadapt_cmd = (
            f"cutadapt "
            f"--length {target_length} "
            f"--output {snakemake.output.fastq} "
            f"{concatenated_fastq}"
        )
        f = open(log_filename, 'at')
        f.write("## COMMAND: "+cutadapt_cmd+"\n")
        f.close()
        shell(cutadapt_cmd)
        
        # Clean up temporary concatenated file if created
        if len(in_fastq_list) > 1:
            command = "rm " + concatenated_fastq
            f = open(log_filename, 'at')
            f.write("## COMMAND: "+command+"\n")
            f.close()
            shell(command)