#############################################################
# wrapper for rule: create_samplesheet
#############################################################
import re
import csv
import pandas as pd



sample_tab = snakemake.params.sample_tab
if snakemake.params.config["run_lane_splitting"] != None:
    column_name = snakemake.params.lane
    column_name = column_name.replace("L0","lane")
    sample_tab = sample_tab[sample_tab[column_name] == True]
sample_tab = sample_tab[sample_tab['demux_setting'] == snakemake.params.demux]

def get_sequenced_barcode_lengths(file_path):
    bioinfo_df = pd.read_csv(file_path)

    # Initialize variables to hold the sequence lengths
    i7_sequence_length = None
    i5_sequence_length = None

    # Loop through the DataFrame to find the relevant information
    for index, row in bioinfo_df.iterrows():
        if 'Barcode' in row.values:
            i7_sequence_length = int(row.values[1])  # Assuming the length is in the second column
        elif 'Dual Barcode' in row.values:
            i5_sequence_length = int(row.values[1])  # Assuming the length is in the second column
        if i7_sequence_length and i5_sequence_length:
            break  # Stop the loop if both lengths are found

    return [i7_sequence_length, i5_sequence_length]

def concat_sequences(i7_seq, i5_seq, i7_len, i5_len):
    if all(char == 'G' for char in i7_seq):
        i7_seq = i7_seq.replace('G', 'N')
    return i5_seq[-i5_len:] + i7_seq[-i7_len:]


sequenced_barcode_lengths = get_sequenced_barcode_lengths(snakemake.input.run_info)



# Apply the function to create a new column with the concatenated sequences
sample_tab['concatenated_sequences'] = sample_tab.apply(lambda row: concat_sequences(row['i7_sequence'],
                                                                                     row['i5_sequence'],
                                                                                     sequenced_barcode_lengths[0],
                                                                                     sequenced_barcode_lengths[1]), axis=1)

# Create a new table with only the sample_name and concatenated_sequences columns
samplesheet_tab = sample_tab[['sample_name_full', 'concatenated_sequences']]

with open(snakemake.output.sample_sheet, 'w') as file:
    samplesheet_tab.to_csv(file, sep='\t', index=False, header=False)

