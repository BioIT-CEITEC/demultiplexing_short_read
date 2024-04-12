######################################
# wrapper for rule: copy_stats
######################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")

out_dir_name = os.path.dirname(snakemake.output)
command = "mkdir -p " + out_dir_name + " && cp demux_info.tsv " + out_dir_name
shell(command)

for filename in snakemake.input.files:
    command = "mkdir -p " + os.path.join(out_dir_name,os.path.dirname(filename))\
              + " && cp -r " + filename + " " + os.path.join(out_dir_name,os.path.dirname(filename))
    shell(command)