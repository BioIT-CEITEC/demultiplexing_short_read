######################################
# wrapper for rule: fastq_prepare_PE
######################################
import os
import re
import gzip
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: fastq_prepare_PE \n##\n")
f.close()

shell.executable("/bin/bash")

command = "mkdir -p " + os.path.dirname(snakemake.output[0])
f = open(log_filename, 'at')
f.write("## CREATE_OUTPUT_DIR: " + command + "\n")
f.close()
shell(command)


def replace_last_occurrence(s, old, new):
    # Reverse the string (so the last occurrence becomes the first)
    reversed_s = s[::-1]

    # Replace the first occurrence of the substring in the reversed string
    reversed_s = reversed_s.replace(old[::-1], new[::-1], 1)

    # Reverse the string back to its original orientation
    s = reversed_s[::-1]

    return s


if os.stat(snakemake.input.in_filename).st_size != 0:
    sample = snakemake.wildcards.sample
    in_filename = snakemake.input.in_filename
    in_filename_R2 = replace_last_occurrence(in_filename, "_R1", "_R2")
    umi = snakemake.params.umi
    run_name = snakemake.params.run_name

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
        f = open(log_filename, 'at')
        f.write("## UMI COMMAND: "+command+"\n")
        f.close()
        shell(command)

        #copy R1
        if "externally_sequenced_fake" in run_name:
            command = "cp -T "+in_filename+" "+ snakemake.output.R1
        else:
            command = "mv -T "+in_filename+" "+ snakemake.output.R1

        f = open(log_filename, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)

        #copy R2
        if "externally_sequenced_fake" in run_name:
            command = "cp -T "+in_filename_R2+" "+ snakemake.output.R2
        else:
            command = "mv -T "+in_filename_R2+" "+ snakemake.output.R2

        f = open(log_filename, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)

    elif umi == "Qiaseq":
        out_R1 = snakemake.output.R1[:-3]
        out_R2 = snakemake.output.R2[:-3]

        command = "(paste <(zcat "+in_filename+") <(zcat "+in_filename_R2+") |"+\
                  " awk '{{ if(NR%4==1) {{split($1,head_R1,\" \"); split($2,head_R2,\" \")}}"+\
                  " else if(NR%4==2) {{umi=substr($2,1,12); print head_R1[1] \"_\" umi \" \" head_R1[2] \"\\n\" $1 > out1;"+\
                  " print head_R2[1] \"_\" umi \" \" head_R2[2] \"\\n\" substr($2,24) > out2}} else if(NR%4==0) {{print $1 > out1;"+\
                  " print substr($2,24) > out2}} else {{print $1 > out1; print $2 > out2}} }}' FS='\\t' out1="+out_R1+" out2="+out_R2+\
                  " && gzip -f "+out_R1+" "+out_R2+") 2>> "+log_filename
        with open(log_filename, 'at') as f:
            f.write("## COMMAND: "+command+"\n")
        shell(command)

        #if "externally_sequenced_fake" in run_name:
        #    command = "cp -T "+in_filename+" "+ snakemake.output.R1
        #else:
        #    command = "mv -T "+in_filename+" "+ snakemake.output.R1

        #f = open(log_filename, 'at')
        #f.write("## COMMAND: "+command+"\n")
        #f.close()
        #shell(command)

        #umi_file = os.path.dirname(snakemake.output.R2)+ "/"+sample+".UMI.fastq"
        #shell("gunzip "+in_filename_R2)
        #in_filename_R2 = in_filename_R2.replace(".gz","")

        #f_R2 = open(snakemake.output.R2.replace(".gz",""), 'w')
        #f_umi = open(umi_file, 'w')

        #with(open(in_filename_R2,'r')) as f_in:
        #    for line_num,line in enumerate(f_in):
        #        if line_num % 2 == 1:
        #            f_umi.write(line[:12]+"\n")
        #            f_R2.write(line[24:])
        #        else:
        #            f_umi.write(line)
        #            f_R2.write(line)

        #f_umi.close()
        #f_R2.close()
        #shell("gzip " + snakemake.output.R2.replace(".gz",""))
        #shell("rm "+in_filename_R2)

    elif umi == "CS_UMI_sep_file":

        umi_file = os.path.dirname(snakemake.output.R2)+ "/"+sample+".UMI.fastq"
        out_R1 = snakemake.output.R1[:-3]
        out_R2 = snakemake.output.R2[:-3]
        
        command = "(paste -d '' <(zcat "+in_filename+" | awk '{{ if(NR%4==2 || NR%4==0) {{print substr($0,1,3); print substr($0,7) > out}} else {{print $0; print $0 > out}} }}' out="+out_R1+") <(zcat "+in_filename_R2+" | awk '{{ if(NR%4==2 || NR%4==0) {{print substr($0,1,3); print substr($0,7) > out}} else {{print \"\"; print $0 > out}} }}' out="+out_R2+") > "+umi_file+" && gzip -f "+out_R1+" "+out_R2+") 2>> "+log_filename
        with open(log_filename, 'at') as f:
            f.write("## COMMAND: "+command+"\n")
        shell(command)


    elif umi == "CS_UMI":
        out_R1 = snakemake.output.R1[:-3]
        out_R2 = snakemake.output.R2[:-3]
        
        command = "(paste <(zcat "+in_filename+") <(zcat "+in_filename_R2+") |"+\
                  " awk '{{ if(NR%4==1) {{split($1,head_R1,\" \"); split($2,head_R2,\" \")}}"+\
                  " else if(NR%4==2) {{umi=substr($1,1,3)substr($2,1,3); print head_R1[1] \"_\" umi \" \" head_R1[2] \"\\n\" substr($1,7) > out1;"+\
                  " print head_R2[1] \"_\" umi \" \" head_R2[2] \"\\n\" substr($2,7) > out2}} else if(NR%4==0) {{print substr($1,7) > out1;"+\
                  " print substr($2,7) > out2}} else {{print $1 > out1; print $2 > out2}} }}' FS='\\t' out1="+out_R1+" out2="+out_R2+\
                  " && gzip -f "+out_R1+" "+out_R2+") 2>> "+log_filename
        with open(log_filename, 'at') as f:
            f.write("## COMMAND: "+command+"\n")
        shell(command)

    elif umi == "TruSight_Oncology":
        out_R1 = snakemake.output.R1[:-3]
        out_R2 = snakemake.output.R2[:-3]

        command = "(paste <(zcat " + in_filename + ") <(zcat " + in_filename_R2 + ") |" + \
                  " awk '{{ if(NR%4==1) {{split($1,head_R1,\" \"); split($2,head_R2,\" \")}}" + \
                  " else if(NR%4==2) {{umi=substr($1,1,6)substr($2,1,6); print head_R1[1] \"_\" umi \" \" head_R1[2] \"\\n\" substr($1,10) > out1;" + \
                  " print head_R2[1] \"_\" umi \" \" head_R2[2] \"\\n\" substr($2,10) > out2}} else if(NR%4==0) {{print substr($1,10) > out1;" + \
                  " print substr($2,10) > out2}} else {{print $1 > out1; print $2 > out2}} }}' FS='\\t' out1=" + out_R1 + " out2=" + out_R2 + \
                  " && gzip -f " + out_R1 + " " + out_R2 + ") 2>> " + log_filename
        with open(log_filename, 'at') as f:
            f.write("## COMMAND: " + command + "\n")
        shell(command)

    else:
        command = "mv -T " + in_filename + " " + snakemake.output.R1
        f = open(log_filename, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)

        command = "mv -T " + in_filename_R2 + " " + snakemake.output.R2
        f = open(log_filename, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)
else:
    command = "touch " + snakemake.output.R1
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    command = "touch " + snakemake.output.R2
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)
