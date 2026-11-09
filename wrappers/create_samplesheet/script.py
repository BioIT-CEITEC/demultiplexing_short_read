#############################################################
# wrapper for rule: create_samplesheet
#############################################################
import re
import xml.etree.ElementTree as ET

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

def read_run_info(run_info_filename):
    """Return [(num_cycles, is_index), ...] for all reads in RunInfo.xml, in read order."""
    root = ET.parse(run_info_filename).getroot()
    reads = sorted(root.iter("Read"), key=lambda r: int(r.get("Number")))
    return [(int(r.get("NumCycles")), r.get("IsIndexedRead") == "Y") for r in reads]


def basemask_to_override_cycles(mask, read_cycles):
    """Convert a bcl2fastq --use-bases-mask (e.g. 'Y*,I8n*,Y10,Y*') into a
    BCL Convert OverrideCycles value (e.g. 'Y151;I8N2;Y10;Y151').
    Masks already in OverrideCycles syntax pass through (validated)."""
    mask = mask.strip()
    if mask == "":
        return ""
    if ":" in mask:
        # bcl2fastq per-lane syntax "1:Y*,I8,Y*" - OverrideCycles is global
        lane_masks = {m.split(":", 1)[1] for m in mask.split() if ":" in m}
        if len(lane_masks) != 1:
            raise ValueError("Per-lane base masks with different values are not supported by OverrideCycles: " + mask)
        mask = lane_masks.pop()

    parts = re.split(r"[,;]", mask)
    if len(parts) != len(read_cycles):
        raise ValueError(f"Base mask '{mask}' has {len(parts)} reads but RunInfo.xml has {len(read_cycles)}")

    out_parts = []
    for part, n_cycles in zip(parts, read_cycles):
        tokens = re.findall(r"([YINUyinu])(\*|\d+)?", part)
        if "".join(l + c for l, c in tokens) != part:
            raise ValueError(f"Cannot parse base mask segment '{part}'")
        if sum(c == "*" for _, c in tokens) > 1:
            raise ValueError(f"More than one '*' in base mask segment '{part}'")

        fixed = sum(int(c) if c.isdigit() else 1 for _, c in tokens if c != "*")
        remaining = n_cycles - fixed
        if remaining < 0:
            raise ValueError(f"Base mask segment '{part}' is longer than the read ({n_cycles} cycles)")

        expanded = []
        for letter, count in tokens:
            n = remaining if count == "*" else (int(count) if count else 1)
            if n > 0:
                expanded.append([letter.upper(), n])
        if remaining > 0 and not any(c == "*" for _, c in tokens):
            # bcl2fastq masked trailing unspecified cycles - make it explicit
            expanded.append(["N", remaining])

        merged = []
        for letter, n in expanded:
            if merged and merged[-1][0] == letter:
                merged[-1][1] += n
            else:
                merged.append([letter, n])
        out_parts.append("".join(f"{l}{n}" for l, n in merged))
    return ";".join(out_parts)


# Main code
run_info_filename = str(snakemake.input.run_info)
run_reads = read_run_info(run_info_filename)
run_index_lengths = [cycles for cycles, is_index in run_reads if is_index]

sample_tab = snakemake.params.sample_tab
sample_tab = sample_tab[sample_tab['demux_setting'] == snakemake.params.demux_setting]

# BCL Convert demultiplexing settings (bcl2fastq took these on the CLI)
demux_setting = sample_tab.iloc[0].to_dict()
basemask = demux_setting.get("illumina_basemask", "")
if not isinstance(basemask, str):  # missing column / NaN
    basemask = ""
override_cycles = basemask_to_override_cycles(basemask, [cycles for cycles, _ in run_reads])

# BCL Convert requires sample sheet index lengths to match the I cycles of
# OverrideCycles (e.g. I8N2 -> 8 bp; an index read turned into Y -> no index)
if override_cycles != "":
    segments = override_cycles.split(";")
    run_index_lengths = [sum(int(n) for letter, n in re.findall(r"([YINU])(\d+)", segments[i]) if letter == "I")
                         for i, (_, is_index) in enumerate(run_reads) if is_index]
# a missing index read in the run means no index sequence in the sample sheet
if len(run_index_lengths) > 0:
    run_index_lengths += [0] * (2 - len(run_index_lengths))

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
    if override_cycles != "":
        writer.write(f'OverrideCycles,{override_cycles}\n')
    writer.write(f'BarcodeMismatchesIndex1,{demux_setting["barcode_mismatches"]}\n')
    if "i5_sequence" in sample_tab and len(run_index_lengths) > 1 and run_index_lengths[1] > 0:
        writer.write(f'BarcodeMismatchesIndex2,{demux_setting["barcode_mismatches"]}\n')
    # keep short (e.g. UMI) reads unpadded, as with the previous bcl2fastq settings
    writer.write('MinimumTrimmedReadLength,2\n')
    writer.write('MaskShortReads,2\n')
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