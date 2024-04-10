#############################################################
# wrapper for rule: create_samplesheet
#############################################################
import re
import csv
import pandas as pd

config = snakemake.params.config
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


#process sample tab
sample_tab = snakemake.params.sample_tab
#add lane info from config
#create lane info table
lane_info = []
# Iterate through each project in the 'libraries' key
for lib_name, lib in config.get('libraries', {}).items():
    if isinstance(lib, dict):
        # Iterate through each key in the individual dictionary
        if config["run_lane_splitting"] != None:
            if lib["lib_lane_splitting"]["lane1"] == True and lib["lib_lane_splitting"]["lane2"] == False:
                lane = "1"
            if lib["lib_lane_splitting"]["lane2"] == True and lib["lib_lane_splitting"]["lane1"] == False:
                lane = "2"
        else:
            lane = "1+2"

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

    if basemask_tab:
        # Write the text
        file.write("[SETTINGS],,\n")
        # Convert the list of dictionaries to a DataFrame and write as CSV
        basemask_tab = pd.DataFrame(basemask_tab)
        basemask_tab.to_csv(file, index=False)

    # Write the second block of text
    file.write("\n[SAMPLES],,\n")
    # Write the pandas DataFrame as CSV
    to_print_sample_tab.to_csv(file, index=False)


