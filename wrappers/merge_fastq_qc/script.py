######################################
# wrapper for rule: merge_fastq_qc
######################################
import subprocess
from os.path import dirname
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: merge_fastq_qc \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

search_path = " ".join([dirname(fastqc_html) for fastqc_html in snakemake.input.html])
if hasattr(snakemake.input, 'minion'):
    search_path += " "+dirname(snakemake.input.minion)
if hasattr(snakemake.input, 'biobloom'):
    search_path += " "+" ".join([dirname(biobloom_tsv) for biobloom_tsv in snakemake.input.biobloom])
if hasattr(snakemake.input, 'sp_det'):
    search_path += " "+dirname(snakemake.input.sp_det)

command = "multiqc -f -n " + snakemake.output.html + " " + search_path + \
              " --cl_config \"{{read_count_multiplier: 0.001, read_count_prefix: 'K', read_count_desc: 'thousands' }}\" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
