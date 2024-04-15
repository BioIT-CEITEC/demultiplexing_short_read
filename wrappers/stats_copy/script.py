######################################
# wrapper for rule: copy_stats
######################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
import xml.etree.ElementTree as ET
import pandas as pd
import json

out_dir_name = os.path.dirname(snakemake.output.nread_json)
sample_tab = snakemake.params.sample_tab
# command = "mkdir -p " + out_dir_name + " && cp demux_info.tsv " + out_dir_name
# shell(command)
if snakemake.params.sequencer_type == "AVITI":

    # Load the CSV file
    index_assignment_path = str(snakemake.params.nread_file[0])
    index_assignment_df = pd.read_csv(index_assignment_path)

    # Summing up the NumPoloniesAssigned for each SampleName
    barcode_counts = index_assignment_df.groupby('SampleName')['NumPoloniesAssigned'].sum().to_dict()

elif snakemake.params.sequencer_type == "MGI":
    file_paths = snakemake.params.nread_file  # Example with one file, add more file names if provided

    # Dictionary to store barcode counts from each file
    barcode_counts = {}

    # Process each file
    for file_path in file_paths:
        # Read the current file
        df = pd.read_csv(file_path, delimiter='\t')

        # Update counts in the dictionary
        for index, row in df.iterrows():
            barcode = row['#Barcode']
            total_count = row['Total']

            if barcode in barcode_counts:
                barcode_counts[barcode] += total_count
            else:
                barcode_counts[barcode] = total_count

else:
    infile = str(snakemake.params.nread_file[0])
    # Load and parse the XML file
    tree = ET.parse(infile)
    root = tree.getroot()

    # Initialize a dictionary to store the sum of BarcodeCounts for each sample
    barcode_counts = {}

    # Loop through XML to extract data
    for sample in root.findall(".//Sample"):
        sample_name = sample.get('name')
        for barcode in sample.findall(".//Barcode"):
            barcode_name = barcode.get('name')
            if barcode_name == "all":
                for lane in barcode.findall("Lane"):
                    count = int(lane.find('BarcodeCount').text)
                    if sample_name in barcode_counts:
                        barcode_counts[sample_name] += count
                    else:
                        barcode_counts[sample_name] = count

# Map summed BarcodeCounts to the DataFrame sample_ID
sample_to_barcode_count = {
    row['sample_ID']: barcode_counts.get(row['sample_name_full'], 0)
    for index, row in sample_tab.iterrows()
}

with open(snakemake.output.nread_json, 'w') as file:
    json.dump(sample_to_barcode_count, file)

for filename in snakemake.params.stats_files:
    command = "mkdir -p " + os.path.join(out_dir_name,os.path.dirname(filename))\
              + " && cp -r " + filename + " " + os.path.join(out_dir_name,os.path.dirname(filename))
    shell(command)