######################################
# wrapper for rule: bcl2fastq
######################################
import os
from snakemake.shell import shell

stats_dir = os.path.dirname(snakemake.input.stats)
out_dir_name = os.path.dirname(snakemake.output.stats)

print("cp -r " + stats_dir + "/ " + out_dir_name)

shell("mkdir -p " + out_dir_name)
shell("cp " + stats_dir + "/* " + out_dir_name)
shell("cp -r " + stats_dir + "/../Reports " + out_dir_name)
