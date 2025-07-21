#############################################################
# wrapper for rule: create_samplesheet
#############################################################
import re
import csv

def write_sample_info(writer, run_index_lengths, sample_tab, lane):
    fake_index_list = ["AAAAAAAAAAAA", "TAAAAAAAAAAA", "TTAAAAAAAAAA", "TTTAAAAAAAAA", "TTTTAAAAAAAA", "TTTTTAAAAAAA", "TTTTTTAAAAAAA"]
    existing_index_check_dict = {}
    dup_index_count = 0

    if lane == 0:
        lane_text = ""
    else:
        lane_text = str(lane) + ","

    if "i5_sequence" in sample_tab:
        if len(run_index_lengths) < 2:
            run_index_lengths = [100, 100]
        for index, row in sample_tab.iterrows():
            i7_seq = row["i7_sequence"][0:run_index_lengths[0]]
            i5_seq = row["i5_sequence"][0:run_index_lengths[1]]
            if i7_seq + i5_seq not in existing_index_check_dict:
                writer.write(lane_text + row["sample_name_full"] + ",,,," +
                             row["i7_name"] + "," + i7_seq + "," +
                             row["i5_name"] + "," + i5_seq + ",,\n")
                existing_index_check_dict[i7_seq + i5_seq] = True
            else:
                writer.write(lane_text + row["sample_name_full"] + ",,,," +
                             row["i7_name"] + "," +
                             fake_index_list[dup_index_count][0:run_index_lengths[0]] + "," +
                             row["i5_name"] + "," +
                             fake_index_list[dup_index_count][0:run_index_lengths[1]] + ",,\n")
                dup_index_count += 1
    else:
        if len(run_index_lengths) < 1:
            run_index_lengths = [100]
        for index, row in sample_tab.iterrows():
            i7_seq = row["i7_sequence"][0:run_index_lengths[0]]
            if i7_seq not in existing_index_check_dict:
                writer.write(lane_text + row["sample_name_full"] + ",,,," +
                             row["i7_name"] + "," + i7_seq + ",,,,\n")
                existing_index_check_dict[i7_seq] = True
            else:
                writer.write(lane_text + row["sample_name_full"] + ",,,," +
                             row["i7_name"] + "," +
                             fake_index_list[dup_index_count][0:run_index_lengths[0]] + ",,,,\n")
                dup_index_count += 1

# Main code
run_info_filename = str(snakemake.input.run_info)

with open(run_info_filename) as f:
    run_index_lengths = re.findall(r"<Read.*NumCycles=\"([0-9]+)\".*IsIndexedRead=\"Y\" />", f.read())
run_index_lengths = [int(index) for index in run_index_lengths]

sample_tab = snakemake.params.sample_tab
sample_tab = sample_tab[sample_tab['demux_setting'] == snakemake.params.demux_setting]

date = str(snakemake.params.date)
with open(snakemake.output.samplesheet_csv, mode='w') as samplesheet_file:
    writer = samplesheet_file  # Use the file object directly for manual writing

    writer.write('[Header]\n')
    writer.write('IEMFileVersion,4\n')
    writer.write(f'Experiment name,{snakemake.params.run_name}\n')
    writer.write(f'Date,{date[2:4]}/{date[4:6]}/20{date[0:2]}\n')
    writer.write('Workflow,GenerateFASTQ\n')
    writer.write('Application,FASTQ Only\n')
    writer.write('Description\n')

    if "i5_sequence" in sample_tab:
        writer.write('Chemistry,Amplicon\n')
    else:
        writer.write('Chemistry,Default\n')

    writer.write('[Reads]\n')
    writer.write(f'{snakemake.params.run_forward_read_length}\n')
    writer.write(f'{snakemake.params.run_reverse_read_length}\n')
    writer.write('[Settings]\n')
    writer.write('[Data]\n')
    if snakemake.params.run_lane_splitting is not None:
        writer.write('Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description\n')

        lane_columns = [col for col in sample_tab.columns if col.startswith('lane')]
        for col in lane_columns:
            lane_number = int(col[4:])
            df_slice = sample_tab[sample_tab[col] == True]
            write_sample_info(writer, run_index_lengths, df_slice, lane_number)
    else:
        writer.write('Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description\n')
        write_sample_info(writer, run_index_lengths, sample_tab, 0)


# In this version, the `writer` object is used directly to write to the file, avoiding the `csv.writer` which introduces unwanted quotes and escape characters. This should produce the desired CSV format without any additional characters.