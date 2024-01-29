######################################
# wrapper for rule: copy_stats
######################################
import os
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: copy_stats \n##\n")
f.close()

#version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
#f = open(snakemake.log.run, 'at')
#f.write("## CONDA: "+version+"\n")
#f.close()

out_dir_name = os.path.dirname(snakemake.output.stats)


command = "mkdir -p " + out_dir_name + " 2>> "+snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

for filename in snakemake.input.files:
    command = "cp " + filename + " " + out_dir_name + "/ 2>> "+snakemake.log.run
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

command = "cp " + snakemake.params.stats_json_file + " " + snakemake.output.stats + " 2>> "+snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)