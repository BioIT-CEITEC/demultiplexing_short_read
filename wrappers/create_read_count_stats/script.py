#############################################################
# wrapper for rule: create_read_count_stats
#############################################################
import subprocess
import json

if snakemake.params.config["library_output"][snakemake.params.library_output]["lib_reverse_read_length"] > 0:
    read_pair_tags = "_R1"
else:
    read_pair_tags = ""

sample_dict = dict(snakemake.params.config["library_output"][snakemake.params.library_output]["samples"])
sample_list, lines_list, dict_list = [], [], []

# prepare input file for counting - only SE or R1 samples
for name in sample_dict.keys():
    tmp = str(snakemake.params.lib_name)+sample_dict[name]["sample_name"]+read_pair_tags+".fastq.gz"
    sample_list.append(tmp)

# count number of lines
for i in range(0,len(sample_list)):
    n_lines = str(subprocess.Popen("gunzip -c " + str(sample_list[i]) + " | wc -l", shell=True,stdout=subprocess.PIPE).communicate()[0], 'utf-8')
    if int(n_lines) != 0:
        n_lines = int(n_lines) / 4
    lines_list.append(int(n_lines))

# create dictionary in json format
for i, name in enumerate(sample_dict.keys()):
    if not "i5_sequence" in sample_dict[name]:
        tmp = {"SampleId": sample_dict[name]["sample_name"], "IndexMetrics": [{"IndexSequence": sample_dict[name]["i7_sequence"]}],"NumberReads": lines_list[i]}
        dict_list.append(tmp)
    else:
        tmp = {"SampleId": sample_dict[name]["sample_name"], "IndexMetrics": [{"IndexSequence": sample_dict[name]["i7_sequence"] + "+" + sample_dict[name]["i5_sequence"]}],"NumberReads": lines_list[i]}
        dict_list.append(tmp)
dictionary = {"ConversionResults": [{"DemuxResults": dict_list}]}


with open(str(snakemake.output), 'w') as f:
    json.dump(dictionary,f,indent = 4)
