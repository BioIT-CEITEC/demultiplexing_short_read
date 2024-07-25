#############################################################
# wrapper for rule: create_samplesheet
#############################################################
import re
import csv

def write_sample_info(writer,run_index_lengths,sample_tab,lane):
    fake_index_list = ["AAAAAAAAAAAA","TAAAAAAAAAAA","TTAAAAAAAAAA","TTTAAAAAAAAA","TTTTAAAAAAAA","TTTTTAAAAAAA","TTTTTTAAAAAAA"]
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
            if not row["i7_sequence"][0:run_index_lengths[0]] + row["i5_sequence"][0:run_index_lengths[1]] in existing_index_check_dict.keys():
                writer.writerow([lane_text + row["sample_name_full"], "", "", "",
                                row["i7_name"],
                                row["i7_sequence"][0:run_index_lengths[0]],
                                row["i5_name"],
                                row["i5_sequence"][0:run_index_lengths[1]], "", ""])
                existing_index_check_dict[row["i7_sequence"][0:run_index_lengths[0]] + row["i5_sequence"][0:run_index_lengths[1]]] = True
            else:
                writer.writerow([lane_text + row["sample_name_full"], "", "", "",
                                 row["i7_name"],
                                 fake_index_list[dup_index_count][0:run_index_lengths[0]],
                                 row["i5_name"],
                                 fake_index_list[dup_index_count][0:run_index_lengths[1]], "", ""])
                dup_index_count = dup_index_count + 1
    else:
        if len(run_index_lengths) < 1:
            run_index_lengths = [100]
        for index, row in sample_tab.iterrows():
            if not row["i7_sequence"][0:run_index_lengths[0]] in existing_index_check_dict.keys():
                writer.writerow([lane_text + row["sample_name_full"], "", "", "",
                                row["i7_name"],
                                row["i7_sequence"][0:run_index_lengths[0]],
                                "",
                                "", "", ""])
                existing_index_check_dict[row["i7_sequence"][0:run_index_lengths[0]]] = True
            else:
                writer.writerow([lane_text + row["sample_name_full"], "", "", "",
                                row["i7_name"],
                                fake_index_list[dup_index_count][0:run_index_lengths[0]],
                                "",
                                "", "", ""])
                dup_index_count = dup_index_count + 1




run_info_filename = str(snakemake.input.run_info)

with open(run_info_filename) as f:
   run_index_lengths = re.findall("<Read.*NumCycles=\"([0-9]+)\".*IsIndexedRead=\"Y\" />", f.read())
run_index_lengths = [int(index) for index in run_index_lengths]
# print(run_index_lengths)

sample_tab = snakemake.params.sample_tab
sample_tab = sample_tab[sample_tab['demux_setting'] == snakemake.params.demux_setting]

date = str(snakemake.params.date)
with open(snakemake.output.samplesheet_csv, mode='w') as samplesheet_file:
    writer = csv.writer(samplesheet_file,quoting=csv.QUOTE_NONE)
    writer.writerow(['[Header]'])
    writer.writerow(['IEMFileVersion', '4'])
    writer.writerow(['Experiment name', snakemake.params.run_name])
    writer.writerow(['Date', (date[2:3] + "/" + date[4:5] + "/20" + date[0:1])])
    writer.writerow(['Workflow', 'GenerateFASTQ'])
    writer.writerow(['Application', 'FASTQ Only'])
    writer.writerow(['Description'])

    if "i5_sequence" in sample_tab:
        writer.writerow(['Chemistry', 'Amplicon'])
    else:
        writer.writerow(['Chemistry', 'Default'])

    writer.writerow(['[Reads]'])
    writer.writerow([str(snakemake.params.run_forward_read_length)])
    writer.writerow([str(snakemake.params.run_reverse_read_length)])
    writer.writerow(['[Settings]'])
    writer.writerow(['[Data]'])
    if snakemake.params.run_lane_splitting != None:
        writer.writerow(
        ['Lane', 'Sample_ID', 'Sample_Name', 'Sample_Plate', 'Sample_Well', 'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',\
         'Sample_Project', 'Description'])

        lane_columns = [col for col in sample_tab.columns if col.startswith('lane')]
        for col in lane_columns:

            lane_number = int(col[4:])
            # Creating a DataFrame slice where the current column's value is True
            df_slice = sample_tab[sample_tab[col] == True]

            # Calling write_sample_info with the column name and the slice
            write_sample_info(writer,run_index_lengths,df_slice,lane_number)

    else:
        writer.writerow(
            ['Sample_ID', 'Sample_Name', 'Sample_Plate', 'Sample_Well', 'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',\
                'Sample_Project', 'Description'])
        write_sample_info(writer,run_index_lengths, sample_tab, 0)
