#############################################################
# wrapper for rule: fastq_mv
#############################################################
import os
from snakemake.shell import shell

shell.executable("/bin/bash")

#create output dir
shell("mkdir -p " + os.path.dirname(snakemake.output.fastq))

in_fastq_list = []
for in_file in snakemake.params.fastq:
    if os.path.isfile(in_file):
        in_fastq_list.append(in_file)

if len(in_fastq_list) == 0:
    shell("touch " + snakemake.output.fastq)
elif len(in_fastq_list) == 1:
    shell("mv " + in_fastq_list[0] + " " + snakemake.output.fastq)
else:
    shell("cat " + " ".join(in_fastq_list) + " > " + snakemake.output.fastq)
