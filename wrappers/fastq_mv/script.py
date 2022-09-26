#############################################################
# wrapper for rule: fastq_mv
#############################################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")

for in_file, out_file in zip(snakemake.params.fastqs_in, snakemake.output.fastqs_out):
    if os.path.isfile(in_file):
        shell("mv " + in_file + " " + out_file)
    else:
        shell("touch " + out_file)
    if os.path.isfile(in_file.replace("_R1_001.fastq.gz","_R2_001.fastq.gz")):
        shell("mv " + in_file.replace("_R1_001.fastq.gz","_R2_001.fastq.gz") + " " + out_file.replace("_R1_001.fastq.gz","_R2_001.fastq.gz"))
    if os.path.isfile(in_file.replace("_R1_001.fastq.gz","_R3_001.fastq.gz")):
        shell("mv " + in_file.replace("_R1_001.fastq.gz","_R3_001.fastq.gz") + " " + out_file.replace("_R1_001.fastq.gz","_R3_001.fastq.gz"))
    if os.path.isfile(in_file.replace("_R1_001.fastq.gz","_R4_001.fastq.gz")):
        shell("mv " + in_file.replace("_R1_001.fastq.gz","_R4_001.fastq.gz") + " " + out_file.replace("_R1_001.fastq.gz","_R4_001.fastq.gz"))
