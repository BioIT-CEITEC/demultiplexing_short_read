#############################################################
# wrapper for rule: create_samplesheet
#############################################################
import re
import csv


run_info_filename = str(snakemake.input.run_info)

with open(run_info_filename) as f:
   run_index_lengths = re.findall("<Read.*NumCycles=\"([0-9]+)\".*IsIndexedRead=\"Y\" />", f.read())
run_index_lengths = [int(index) for index in run_index_lengths]
print(run_index_lengths)

sample_tab = snakemake.params.sample_tab


date = str(snakemake.params.date)
with open(snakemake.output.samplesheet_csv, mode='w') as samplesheet_file:
    writer = csv.writer(samplesheet_file)
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
    writer.writerow(
        ['Sample_ID', 'Sample_Name', 'Sample_Plate', 'Sample_Well', 'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
         'Sample_Project', 'Description'])

    fake_index_list = ["AAAAAAAAAAAA","TAAAAAAAAAAA","TTAAAAAAAAAA","TTTAAAAAAAAA","TTTTAAAAAAAA","TTTTTAAAAAAA","TTTTTTAAAAAAA"]
    existing_index_check_dict = {}
    dup_index_count = 0

    if "i5_sequence" in sample_tab:
        if len(run_index_lengths) < 2:
            run_index_lengths = [100, 100]
        for index, row in sample_tab.iterrows():
            if not row["i7_sequence"][0:run_index_lengths[0]] + row["i5_sequence"][0:run_index_lengths[1]] in existing_index_check_dict.keys():
                writer.writerow([row["sample_name"], "", "", "",
                                row["i7_name"],
                                row["i7_sequence"][0:run_index_lengths[0]],
                                row["i5_name"],
                                row["i5_sequence"][0:run_index_lengths[1]], "", ""])
                existing_index_check_dict[row["i7_sequence"][0:run_index_lengths[0]] + row["i5_sequence"][0:run_index_lengths[1]]] = True
            else:
                writer.writerow([row["sample_name"], "", "", "",
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
                writer.writerow([row["sample_name"], "", "", "",
                                row["i7_name"],
                                row["i7_sequence"][0:run_index_lengths[0]],
                                "",
                                "", "", ""])
                existing_index_check_dict[row["i7_sequence"][0:run_index_lengths[0]]] = True
            else:
                writer.writerow([row["sample_name"], "", "", "",
                                row["i7_name"],
                                fake_index_list[dup_index_count][0:run_index_lengths[0]],
                                "",
                                "", "", ""])
                dup_index_count = dup_index_count + 1

