######################################
# wrapper for rule: fastq_prepare_PE
######################################
import os
import re
import sys
import math
import subprocess
import gzip
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: fastq_prepare_PE \n##\n")
f.close()

shell.executable("/bin/bash")



command = "mkdir -p " + os.path.dirname(snakemake.output[0])
f = open(log_filename, 'at')
f.write("## CREATE_OUTPUT_DIR: "+command+"\n")
f.close()
shell(command)

sample = snakemake.wildcards.sample
in_filename = snakemake.input.in_filename
in_filename_R2 = in_filename.replace("_R1","_R2")
umi = snakemake.params.umi
run_name = snakemake.params.run_name

if umi == "CORALL":
    command = "umi_tools extract --extract-method=string"+\
            " --bc-pattern=NNNNNNNNNNNN" +\
            " --stdin=" + in_filename +\
            " --stdout=" + snakemake.output.fastq +\
            " -L " + snakemake.log.run

    f = open(snakemake.log.run, 'at')
    f.write("## UMI COMMAND: "+command+"\n")
    f.close()
    shell(command)

    shell("rm " + in_filename)

elif umi == "IDT":
    out_R1 = gzip.open(snakemake.output.R1, 'wt')
    out_R2 = gzip.open(snakemake.output.R2, 'wt')

    in_filename_UMI = in_filename_R2
    in_filename_R2 = re.sub("_R2_","_R3_",in_filename_R2)

    with gzip.open(in_filename,'rt') as R1, gzip.open(in_filename_R2,'rt') as R2, gzip.open(in_filename_UMI,'rt') as UMI:
        i = 0
        for R1_line, R2_line, UMI_line in zip(R1, R2, UMI):
            i += 1
            if i % 4 == 1:
                header_R1 = R1_line.strip()
                header_R2 = R2_line.strip()
            elif i % 4 == 2:
                out_R1.write(header_R1.split(" ")[0] + "_" + UMI_line.strip() + " " + header_R1.split(" ")[1] + "\n" + R1_line)
                out_R2.write(header_R2.split(" ")[0] + "_" + UMI_line.strip() + " " + header_R2.split(" ")[1] + "\n" + R2_line)
            elif i % 4 == 0:
                out_R1.write(R1_line)
                out_R2.write(R2_line)
            else:
                out_R1.write(R1_line)
                out_R2.write(R2_line)

elif umi == "BRONCO":
    out_R1 = gzip.open(snakemake.output.R1, 'wt')
    out_R2 = gzip.open(snakemake.output.R2, 'wt')

    in_filename_UMI = in_filename_R2
    in_filename_R2 = re.sub("_R2_","_R3_",in_filename_R2)

    with gzip.open(in_filename,'rt') as R1, gzip.open(in_filename_R2,'rt') as R2, gzip.open(in_filename_UMI,'rt') as UMI:
        i = 0
        for R1_line, R2_line, UMI_line in zip(R1, R2, UMI):
            i += 1
            if i % 4 == 1:
                header_R1 = R1_line.strip()
                header_R2 = R2_line.strip()
            elif i % 4 == 2:
                out_R1.write(header_R1.split(" ")[0] + "_" + UMI_line.strip() + " " + header_R1.split(" ")[1] + "\n" + R1_line)
                out_R2.write(header_R2.split(" ")[0] + "_" + UMI_line.strip() + " " + header_R2.split(" ")[1] + "\n" + R2_line)
            elif i % 4 == 0:
                out_R1.write(R1_line)
                out_R2.write(R2_line)
            else:
                out_R1.write(R1_line)
                out_R2.write(R2_line)

elif umi == "LYNX":
    umi_file = os.path.dirname(snakemake.output.R2)+ "/"+sample+".UMI.fastq"
    umi_file_in = in_filename_R2
    in_filename_R2 = re.sub("_R2_","_R3_",in_filename_R2)

    command = "gunzip -c "+umi_file_in+" > "+umi_file
    f = open(snakemake.log.run, 'at')
    f.write("## UMI COMMAND: "+command+"\n")
    f.close()
    shell(command)

    #copy R1
    if "externally_sequenced_fake" in run_name:
        command = "cp -T "+in_filename+" "+ snakemake.output.R1
    else:
        command = "mv -T "+in_filename+" "+ snakemake.output.R1

    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

    #copy R2
    if "externally_sequenced_fake" in run_name:
        command = "cp -T "+in_filename_R2+" "+ snakemake.output.R2
    else:
        command = "mv -T "+in_filename_R2+" "+ snakemake.output.R2

    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

elif umi == "Qiaseq":
    if "externally_sequenced_fake" in run_name:
        command = "cp -T "+in_filename+" "+ snakemake.output.R1
    else:
        command = "mv -T "+in_filename+" "+ snakemake.output.R1

    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

    umi_file = os.path.dirname(snakemake.output.R2)+ "/"+sample+".UMI.fastq"
    shell("gunzip "+in_filename_R2)
    in_filename_R2 = in_filename_R2.replace(".gz","")

    f_R2 = open(snakemake.output.R2.replace(".gz",""), 'w')
    f_umi = open(umi_file, 'w')

    with(open(in_filename_R2,'r')) as f_in:
        for line_num,line in enumerate(f_in):
            if line_num % 2 == 1:
                f_umi.write(line[:12]+"\n")
                f_R2.write(line[24:])
            else:
                f_umi.write(line)
                f_R2.write(line)

    f_umi.close()
    f_R2.close()
    shell("gzip " + snakemake.output.R2.replace(".gz",""))
    shell("rm "+in_filename_R2)

elif umi == "CS_UMI_sep_file":

    umi_file = os.path.dirname(snakemake.output.R2)+ "/"+sample+".UMI.fastq"
    f_umi = open(umi_file, 'w')
    out_R1 = gzip.open(snakemake.output.R1, 'wt')
    out_R2 = gzip.open(snakemake.output.R2, 'wt')

    with gzip.open(in_filename,'rt') as R1, gzip.open(in_filename_R2,'rt') as R2:
        i = 0
        for R1_line, R2_line in zip(R1, R2):
            i += 1
            if i % 4 == 1:
                header_R1 = R1_line.strip()
                header_R2 = R2_line.strip()
            elif i % 4 == 2:
                out_R1.write(header_R1 + "\n" + R1_line[6:])
                out_R2.write(header_R2 + "\n" + R2_line[6:])
                f_umi.write(header_R1 + "\n" + R1_line.strip()[0:3] + R2_line.strip()[0:3] + "\n")
            elif i % 4 == 0:
                out_R1.write(R1_line[6:])
                out_R2.write(R2_line[6:])
                f_umi.write(R1_line.strip()[0:3] + R2_line.strip()[0:3] + "\n")
            else:
                out_R1.write(R1_line)
                out_R2.write(R2_line)
                f_umi.write(R1_line)

elif umi == "CS_UMI":
    out_R1 = gzip.open(snakemake.output.R1, 'wt')
    out_R2 = gzip.open(snakemake.output.R2, 'wt')

    with gzip.open(in_filename,'rt') as R1, gzip.open(in_filename_R2,'rt') as R2:
        i = 0
        for R1_line, R2_line in zip(R1, R2):
            i += 1
            if i % 4 == 1:
                header_R1 = R1_line.strip()
                header_R2 = R2_line.strip()
            elif i % 4 == 2:
                out_R1.write(header_R1.split(" ")[0] + "_" + R1_line.strip()[0:3] + R2_line.strip()[0:3] + " " + header_R1.split(" ")[1] + "\n" + R1_line[6:])
                out_R2.write(header_R2.split(" ")[0] + "_" + R1_line.strip()[0:3] + R2_line.strip()[0:3] + " " + header_R2.split(" ")[1] + "\n" + R2_line[6:])
            elif i % 4 == 0:
                out_R1.write(R1_line[6:])
                out_R2.write(R2_line[6:])
            else:
                out_R1.write(R1_line)
                out_R2.write(R2_line)


else:
    #copy R1
    # if "externally_sequenced_fake" in run_name:
    #     command = "cp -T "+in_filename+" "+ snakemake.output.R1
    # else:
    #     command = "mv -T "+in_filename+" "+ snakemake.output.R1
    command = "mv -T " + in_filename + " " + snakemake.output.R1
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

    #copy R2
    # if "externally_sequenced_fake" in run_name:
    #     command = "cp -T "+in_filename_R2+" "+ snakemake.output.R2
    # else:
    #     command = "mv -T "+in_filename_R2+" "+ snakemake.output.R2
    command = "mv -T " + in_filename_R2 + " " + snakemake.output.R2
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

# else:
#     #merge input fastq files
#     command = " cat " + " ".join(sorted(snakemake.params.rep_fastq_R1)) + " > " + snakemake.output.R1
#     f = open(snakemake.log.run, 'at')
#     f.write("## COMMAND: "+command+"\n")
#     f.close()
#     shell(command)
#     command = " cat " + " ".join(sorted(snakemake.params.rep_fastq_R2)) + " > " + snakemake.output.R2
#     f = open(snakemake.log.run, 'at')
#     f.write("## COMMAND: "+command+"\n")
#     f.close()
#     shell(command)
