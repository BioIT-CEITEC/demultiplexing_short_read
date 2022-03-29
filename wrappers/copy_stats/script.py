######################################
# wrapper for rule: bcl2fastq
######################################
import os
from snakemake.shell import shell

stats_dir = os.path.dirname(snakemake.input)
out_dir_name = os.path.dirname(snakemake.output)

shell("mkdir -p " + out_dir_name)
shell("cp -r " + stats_dir + "/ " + out_dir_name)
shell("cp -r " + stats_dir + "/../Reports " + out_dir_name)
