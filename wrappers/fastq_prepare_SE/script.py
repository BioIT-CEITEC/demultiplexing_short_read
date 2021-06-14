######################################
# wrapper for rule: fastq_prepare_SE
######################################
import os
import re
import sys
import math
import time
from glob import glob
import subprocess
import gzip
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: fastq_prepare_SE \n##\n")
f.close()

shell.executable("/bin/bash")


command = "mkdir -p " + os.path.dirname(snakemake.output[0])
f = open(log_filename, 'at')
f.write("## CREATE_OUTPUT_DIR: "+command+"\n")
f.close()
shell(command)

sample = snakemake.wildcards.sample
in_filename = snakemake.input.in_filename
umi = snakemake.params.umi
run_name = snakemake.params.run_name


f = open(log_filename, 'at')
f.write("## UMI SETTING: "+umi+"\n")
f.close()


if umi == "CORALL":
    command = "umi_tools extract --extract-method=string"+\
            " --bc-pattern=NNNNNNNNNNNN" +\
            " --stdin=" + in_filename +\
            " --stdout=" + snakemake.output.fastq +\
            " -L " + log_filename

    f = open(log_filename, 'at')
    f.write("## UMI COMMAND: "+command+"\n")
    f.close()
    shell(command)

    shell("rm " + in_filename)

elif umi == "custom_umi":

#        command = "umi_tools extract --extract-method=string"+\
#                " --bc-pattern=NNNNNNNN" +\
#                " --stdin=" + in_filename +\
#                " --stdout=" + snakemake.output.fastq +\
#                " -L " + log_filename
#
#        f = open(log_filename, 'at')
#        f.write("## UMI COMMAND: "+command+"\n")
#        f.close()
#        shell(command)


     out_SE = gzip.open(snakemake.output.fastq, 'wt')

     in_filename_R2 = re.sub("_R1_","_R2_",in_filename)

     with gzip.open(in_filename,'rt') as R1, gzip.open(in_filename_R2,'rt') as R2:
         i = 0
         cut_length = 0
         for R1_line, R2_line in zip(R1, R2):
             i += 1
             if i % 4 == 1:
                 header_R1 = R1_line.strip()
             elif i % 4 == 2:
                 # R2_seq = Seq(R2_line.strip())
                 # R2_seq = str(R2_seq.reverse_complement())
                 # R2_seq = R2_line.strip()[::-1]
                 # cut_length = len(re.sub('[AN]+$|[TN]+$', '', R2_seq))
                 # out_SE.write(header_R1.split(" ")[0] + "_" + R1_line.strip()[0:6] + " " + header_R1.split(" ")[1] + "\n" + R1_line.strip()[6:] + R2_seq[:cut_length] + "\n")

                 # out_SE.write(header_R1.split(" ")[0] + "_" + R1_line.strip()[0:3] + R2_line.strip()[0:3] + " " + header_R1.split(" ")[1] + "\n" + R1_line[3:])

                 ## Chip_Tanja
                 out_SE.write(header_R1.split(" ")[0] + "_" + R2_line.strip() + " " + header_R1.split(" ")[1] + "\n" + R1_line)
             # elif i % 4 == 0:
             #     out_SE.write(R1_line[3:])
             else:
                 out_SE.write(R1_line)
elif umi == "CS_UMI":
    out_SE = gzip.open(snakemake.output.fastq, 'wt')

    in_filename_R2 = re.sub("_R1_","_R2_",in_filename)

    with gzip.open(in_filename,'rt') as R1, gzip.open(in_filename_R2,'rt') as R2:
        i = 0
        for R1_line, R2_line in zip(R1, R2):
            i += 1
            if i % 4 == 1:
                header_R1 = R1_line.strip()
            elif i % 4 == 2:
                out_SE.write(header_R1.split(" ")[0] + "_" + R1_line.strip()[0:3] + R2_line.strip()[0:3] + " " + header_R1.split(" ")[1] + "\n" + R1_line[6:])
            elif i % 4 == 0:
                out_SE.write(R1_line[6:])
            else:
                out_SE.write(R1_line)

elif umi == "Quantseq FWD":
    command = "umi_tools extract --extract-method=string"+\
            " --bc-pattern=NNNNNN" +\
            " --stdin=" + in_filename +\
            " --stdout=" + snakemake.output.fastq +\
            " -L " + log_filename

    f = open(log_filename, 'at')
    f.write("## UMI COMMAND: "+command+"\n")
    f.close()
    shell(command)

    # shell("rm " + in_filename)

elif umi == "IDT":
    umi_file_in = re.sub("_R1_","_R2_",in_filename)

    command = "gunzip -c "+umi_file_in+" > "+umi_file
    f = open(log_filename, 'at')
    f.write("## UMI COMMAND: "+command+"\n")
    f.close()
    shell(command)

elif umi == "BRB":
    constant_index_name = re.sub(".*;","",in_filename)
    in_filename = re.sub(".gz;.*","",in_filename)
    library_name = re.sub("^[0-9]+_","",snakemake.params.lib_name)
    temp_dir = os.path.dirname(in_filename)
    print(in_filename)
    print(snakemake.params.is_first)


    if not snakemake.params.is_first:
        while not os.path.exists(os.path.join(temp_dir, constant_index_name + "pool_" + library_name + ".files_ready")):
            time.sleep(1)

    else:
        try:
            command = "gunzip "+ os.path.join(temp_dir, constant_index_name + "pool_" + library_name + "_S*")
            f = open(log_filename, 'at')
            f.write("## COMMAND: "+command+"\n")
            f.close()
            shell(command)

            UMI_file = os.path.join(temp_dir,glob(os.path.join(temp_dir, constant_index_name + "pool_" + library_name + "*R1*"))[0])
            SEQ_file = os.path.join(temp_dir,glob(os.path.join(temp_dir, constant_index_name + "pool_" + library_name + "*R2*"))[0])

            with open(UMI_file) as UMI, open(SEQ_file) as SEQ:
                i = 0
                header=""
                key=""
                outputs_dict = dict()
                for x, y in zip(UMI, SEQ):
                    i += 1
                    if i % 4 == 1:
                        header = y.strip()
                    elif i % 4 == 2:
                        header = header.split(" ")[0] + "_" + x.strip()[6:16] + " " + header.split(" ")[1] + "\n"
                        key = x.strip()[:6]
                        if key in outputs_dict:
                             # append the new number to the existing array at this slot
                             outputs_dict[key].append(header)
                        else:
                             # create a new array in this slot
                             outputs_dict[key] = [header]
                        outputs_dict[key].append(y)
                    else:
                        outputs_dict[key].append(y)

            for key in outputs_dict.keys():
                out_file = os.path.join(temp_dir,key + "_" + library_name + ".fastq")
                with open(out_file,"w") as out:
                    out.write("".join(outputs_dict[key]))
        except:
            pass

        command = "touch "+ os.path.join(temp_dir, constant_index_name + "pool_" + library_name + ".files_ready")
        f = open(log_filename, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)

    if os.path.isfile(in_filename):
        command = "gzip "+ in_filename
        f = open(log_filename, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)

        command = "mv -T "+in_filename+".gz "+ snakemake.output.fastq
        f = open(log_filename, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)
# elif umi == "Qiaseq miRNA":
#     command = "umi_tools extract --extract-method=regex"+\
#             " --bc-pattern='.*(?P<discard_1>AACTGTAGGCACCATCAAT)(?P<umi_1>.{{12}})(?P<discard_2>.*)'" +\
#             " --stdin=" + in_filename +\
#             " --stdout=" + snakemake.output.fastq +\
#             " -L " + log_filename
#
#     f = open(log_filename, 'at')
#     f.write("## UMI COMMAND: "+command+"\n")
#     f.close()
#     shell(command)
#
#     shell("rm " + in_filename)
else:
    #copy R1
    # if "externally_sequenced_fake" in run_name:
    #     command = "cp -T "+in_filename+" "+ snakemake.output.fastq
    # else:
    command = "mv -T "+in_filename+" "+ snakemake.output.fastq

    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
# else:
#     #merge input fastq files
#     command = " cat " + " ".join(sorted(snakemake.params.rep_fastq)) + " > " + snakemake.output.fastq
#     f = open(log_filename, 'at')
#     f.write("## COMMAND: "+command+"\n")
#     f.close()
#     shell(command)
