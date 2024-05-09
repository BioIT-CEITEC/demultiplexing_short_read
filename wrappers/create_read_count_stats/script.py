#############################################################
# wrapper for rule: create_read_count_stats
#############################################################
import subprocess
import json

def count_lines(fastq_filename):
    n_lines = str(subprocess.Popen("gunzip -c " + str(fastq_filename) + " | wc -l", shell=True,stdout=subprocess.PIPE).communicate()[0], 'utf-8')
    if int(n_lines) != 0:
        n_lines = int(n_lines) / 4
    else:
        n_lines = 0
    return(n_lines)

sample_dict = dict(snakemake.params.config["library_output"][str(snakemake.params.lib_name)]["samples"])

sample_to_barcode_count = {}

# prepare input file for counting - only SE or R1 samples
for name in sample_dict.keys():
    fastq_filename = str(snakemake.params.lib_name)+"/raw_fastq/"+sample_dict[name]["sample_name"]+"R1.fastq.gz"
    sample_to_barcode_count[name] = count_lines(fastq_filename)

with open(str(snakemake.output), 'w') as file:
    json.dump(sample_to_barcode_count, file)
