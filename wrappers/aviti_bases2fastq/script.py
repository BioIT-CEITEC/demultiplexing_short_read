######################################
# wrapper for rule: bcl2fastq
######################################
import os
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: bases2fastq \n##\n")
f.close()

#download and install tool
tmp_tool_dir = os.path.join(snakemake.params.tmp_dir,"bases2fastq_exec")
if not os.path.exists(tmp_tool_dir):
    os.makedirs(tmp_tool_dir)

command = "wget https://bases2fastq-release.s3.amazonaws.com/bases2fastq-latest.tar.gz -O " + tmp_tool_dir + "/bases2fastq-latest.tar.gz >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "tar -xvf " + tmp_tool_dir + "/bases2fastq-latest.tar.gz --directory " + tmp_tool_dir + " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "pip3 install numpy==1.* bs4==0.* >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "pip3 install 'bokeh>=2.3,<3' >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

version = str(subprocess.Popen(tmp_tool_dir + "/bases2fastq --version", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## INFO: "+version+"\n")
f.close()

#run the demultiplexing
run_dir = os.path.dirname(snakemake.input.run_complete_check[0])

tmp_run_data = os.path.join(snakemake.params.tmp_dir,"run_tmp_data")
if not os.path.exists(tmp_run_data):
    os.makedirs(tmp_run_data)

command = "rsync -rt " + run_dir + "/* " + tmp_run_data + " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

fastq_output_dir = os.path.dirname(snakemake.input.run_manifest)
if fastq_output_dir == "":
    fastq_output_dir = "."

command = tmp_tool_dir + "/bases2fastq " + tmp_run_data \
                 + " " + fastq_output_dir \
                 + " -r " + snakemake.input.run_manifest \
                 + " -p " + str(snakemake.threads) \
                 + " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# command = "multiqc -f -z -n "+snakemake.output.html+\
#           " "+snakemake.output.stats+\
#           " >> "+log_filename+" 2>&1"
# f = open(log_filename, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)
