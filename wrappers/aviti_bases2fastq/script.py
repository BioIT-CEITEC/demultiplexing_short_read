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

def format_basemask_settings():
    sample_tab = snakemake.params.sample_tab

    #create_basemask_tab
    basemask_tab = []
    # Parameters to extract
    setting_dict = sample_tab.iloc[0].to_dict()

    parameters = ["R1FastQMask", "R2FastQMask", "I1Mask", "I2Mask", "UmiMask"]
    # Iterate through each project in the 'libraries' key

    for setting_name, value in setting_dict.items():
        if setting_name.startswith('AVITI'):
            basemask_tab.append({
                "SettingName": setting_name.replace("AVITI_", ""),
                "Value": str(value)
            })
        if setting_name == "barcode_mismatches":
            basemask_tab.append({
                "SettingName": "I1MismatchThreshold",
                "Value": str(value)
            })
            basemask_tab.append({
                "SettingName": "I2MismatchThreshold",
                "Value": str(value)
            })
    # Filter out entries where the value is not an empty string
    basemask_tab = [entry for entry in basemask_tab if entry["Value"] != ""]

    # Format each entry into the '--settings "Key,Value"' format
    settings_strs = []
    for setting in basemask_tab:
        setting_str = f'--settings "{setting["SettingName"]},{setting["Value"]}"'
        settings_strs.append(setting_str)

    
    settings_strs.append(setting_dict["demultiplex_additional_options"])
    # Join all settings strings into a single command line argument
    command_line_arg = " ".join(settings_strs)

    if len(command_line_arg) == 0:
        command_line_arg = setting_dict["demultiplex_additional_options"]
    elif setting_dict["demultiplex_additional_options"] == "":
        command_line_arg += " " + str(setting_dict["demultiplex_additional_options"])

    return command_line_arg


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
run_dir = snakemake.params.run_dir

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

command_line_arg = format_basemask_settings()

command = tmp_tool_dir + "/bases2fastq " + tmp_run_data \
                 + " " + fastq_output_dir \
                 + " -r " + snakemake.input.run_manifest \
                 + " -p " + str(snakemake.threads) \
                 + " " + command_line_arg \
                 + " >> " + log_filename + " 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

if not os.path.isfile("demux_info.tsv"):
    with open("demux_info.tsv", 'w') as file:
        # Write the specified text to the file
        file.write("demux_id\tlane\trun_command\tlibraries\n")

with open("demux_info.tsv", 'w') as file:
    # Write the specified text to the file
    file.write(snakemake.wildcards.demux_setting+"\tall\t"+command+"\t"+";".join(snakemake.params.sample_tab['library'].unique().tolist())+"\n")

command = "touch " + snakemake.output.demultiplex_complete
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
