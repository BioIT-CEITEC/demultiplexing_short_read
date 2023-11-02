######################################
# wrapper for rule: fastq_prepare_SE
######################################
import os
import re
import time
from glob import glob
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

if os.stat(snakemake.input.in_filename).st_size != 0:
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
                     ## Chip_Tanja
                     out_SE.write(header_R1.split(" ")[0] + "_" + R2_line.strip() + " " + header_R1.split(" ")[1] + "\n" + R1_line)
                 else:
                     out_SE.write(R1_line)
                     
    elif umi == "CS_UMI":
        out_SE = snakemake.output.fastq[:-3]

        in_filename_R2 = re.sub("_R1_","_R2_",in_filename)
        
        command = "(paste <(zcat "+in_filename+") <(zcat "+in_filename_R2+") |"+\
                  " awk '{{ if(NR%4==1) {{split($1,head_R1,\" \")}}"+\
                  " else if(NR%4==2) {{umi=substr($1,1,3)substr($2,1,3); print head_R1[1] \"_\" umi \" \" head_R1[2] \"\\n\" substr($1,7) > out}}"+\
                  " else if(NR%4==0) {{print substr($1,7) > out}} else {{print $1 > out}} }}' FS='\\t' out="+out_SE+\
                  " && gzip -f "+out_SE+") 2>> "+log_filename
        with open(log_filename, 'at') as f:
            f.write("## COMMAND: "+command+"\n")
        shell(command)

    elif umi == "TruSight_Oncology":
        out_SE = snakemake.output.fastq[:-3]

        in_filename_R2 = re.sub("_R1_","_R2_",in_filename)
        
        command = "(paste <(zcat "+in_filename+") <(zcat "+in_filename_R2+") |"+\
                  " awk '{{ if(NR%4==1) {{split($1,head_R1,\" \")}}"+\
                  " else if(NR%4==2) {{umi=substr($1,1,6)substr($2,1,6); print head_R1[1] \"_\" umi \" \" head_R1[2] \"\\n\" substr($1,10) > out}}"+\
                  " else if(NR%4==0) {{print substr($1,10) > out}} else {{print $1 > out}} }}' FS='\\t' out="+out_SE+\
                  " && gzip -f "+out_SE+") 2>> "+log_filename
        with open(log_filename, 'at') as f:
            f.write("## COMMAND: "+command+"\n")
        shell(command)

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

    else:

        command = "mv -T "+in_filename+" "+ snakemake.output.fastq

        f = open(log_filename, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)

else:
    command = "touch " + snakemake.output.fastq
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

