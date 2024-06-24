#############################################################
# wrapper for rule: fastq_mv
#############################################################
import os
from snakemake.shell import shell
from datetime import datetime

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

if snakemake.params.run_sequencer_type == "MGI":


    input_fastq = snakemake.output.fastq
    output_fastq = snakemake.output.fastq + ".tmp"
    date_id = snakemake.params.date_id

    # Construct the bash command
    # command = f"zcat {input_fastq} | awk -v flowcell_ID={date_id} '{{if (NR % 4 == 1) {{$0 = \"@\" flowcell_ID substr($0, 2)}}; print}}' | gzip > {output_fastq}"

    # command = (
    #         "zcat " + input_fastq + " | "
    #         + "sed 's/^@FS/@{}FS/'".format(date_id) + " | "
    #         + "gzip > " + output_fastq
    # )

    command = (
            "zcat " + input_fastq + " | "
            + "awk 'NR % 4 == 1 {sub(/^@/, \"@" + date_id + "\")}; {print}' | "
            + "gzip > " + output_fastq
    )
    shell(command)

    command = "rm " + snakemake.output.fastq
    shell(command)

    command = "mv " + snakemake.output.fastq + ".tmp " + snakemake.output.fastq
    shell(command)
    # shell("mv " + snakemake.output.fastq + ".tmp " + snakemake.output.fastq)