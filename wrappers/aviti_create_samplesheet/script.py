#############################################################
# wrapper for rule: create_samplesheet
#############################################################
import re
import csv
import pandas as pd

config = snakemake.params.config

#create_basemask_tab
basemask_tab = []
# Parameters to extract
parameters = ["R1FastQMask", "R2FastQMask", "I1Mask", "I2Mask", "UmiMask"]
# Iterate through each project in the 'libraries' key
for project, libraries in config.get('libraries', {}).items():
    if isinstance(libraries, dict):
        # Iterate through each key in the individual dictionary
        for setting_name, value in libraries.items():
            if setting_name in parameters:
                basemask_tab.append({
                    "SettingName": setting_name,
                    "Value": str(value),
                    "Project": project
                })
            if setting_name == "barcode_mismatches":
                basemask_tab.append({
                    "SettingName": "I1MismatchThreshold",
                    "Value": str(value),
                    "Project": project
                })
                basemask_tab.append({
                    "SettingName": "I2MismatchThreshold",
                    "Value": str(value),
                    "Project": project
                })
# Filter out entries where the value is not an empty string
basemask_tab = [entry for entry in basemask_tab if entry["Value"] != ""]


#process sample tab
sample_tab = snakemake.params.sample_tab
#add lane info from config
#create lane info table
lane_info = []
# Iterate through each project in the 'libraries' key
for lib_name, lib in config.get('libraries', {}).items():
    if isinstance(lib, dict):
        # Iterate through each key in the individual dictionary
        lane = "1+2"
        if lib["select_lanes"] == True:
            if lib["use_lane1"] == True and lib["use_lane2"] == False:
                lane = "1"
            if lib["use_lane2"] == True and lib["use_lane1"] == False:
                lane = "2"
        lane_info.append({
            "library": lib_name,
            "Lane": lane,
        })
lane_info = pd.DataFrame(lane_info)
sample_tab = pd.merge(sample_tab, lane_info, on='library', how='inner')

if "i5_sequence" in sample_tab:
    to_print_sample_tab = sample_tab[['sample_name', 'i7_sequence','i5_sequence','Lane','library']]
    to_print_sample_tab = to_print_sample_tab.rename(columns={'sample_name': 'SampleName', 'i7_sequence': 'Index1','i5_sequence': 'Index2', 'Lane': 'Lane','library': 'Project'})
else:
    to_print_sample_tab = sample_tab[['sample_name', 'i7_sequence', 'Lane', 'library']]
    to_print_sample_tab = to_print_sample_tab.rename(columns={'sample_name': 'SampleName', 'i7_sequence': 'Index1','Lane': 'Lane', 'library': 'Project'})

#print to csv file
with open(snakemake.output.run_manifest, 'a') as file:

    # if basemask_tab:
    #     # Write the text
    #     file.write("[SETTINGS],,\n")
    #     # Convert the list of dictionaries to a DataFrame and write as CSV
    #     basemask_tab = pd.DataFrame(basemask_tab)
    #     basemask_tab.to_csv(file, index=False)

    # Write the second block of text
    file.write("\n[SAMPLES],,\n")
    # Write the pandas DataFrame as CSV
    to_print_sample_tab.to_csv(file, index=False)

#
#
# with open(snakemake.output.samplesheet_csv, mode='w') as samplesheet_file:
#     writer = csv.writer(samplesheet_file)
#     writer.writerow(['[SETTINGS],,'])
#     writer.writerow(['IEMFileVersion', '4'])
#     writer.writerow(['Experiment name', snakemake.params.run_name])
#     writer.writerow(['Workflow', 'GenerateFASTQ'])
#     writer.writerow(['Application', 'FASTQ Only'])
#     writer.writerow(['Description'])
#
#     if "i5_sequence" in sample_tab:
#         writer.writerow(['Chemistry', 'Amplicon'])
#     else:
#         writer.writerow(['Chemistry', 'Default'])
#
#     writer.writerow(['[Reads]'])
#     writer.writerow([str(snakemake.params.run_forward_read_length)])
#     writer.writerow([str(snakemake.params.run_reverse_read_length)])
#     writer.writerow(['[Settings]'])
#     writer.writerow(['[Data]'])
#     writer.writerow(
#         ['Sample_ID', 'Sample_Name', 'Sample_Plate', 'Sample_Well', 'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
#          'Sample_Project', 'Description'])
#
#
