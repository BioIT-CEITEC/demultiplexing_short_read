# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Snakemake workflow for demultiplexing short-read sequencing data from Illumina, Aviti, and MGI platforms. Converts raw sequencing output (BCL files, Bases, etc.) into sample-specific FASTQ files.

## Key Configuration

The workflow is configured via `config.json` (not in repo - created at runtime). Key parameters:
- `run_dir`: Path to sequencing run output
- `run_sequencer_type`: One of "ILLUMINA", "AVITI", or "MGI"
- `libraries`: Dictionary of library configurations with samples, read lengths, and platform-specific masks
- `run_lane_splitting`: Optional lane splitting configuration

## Workflow Structure

**Main entry point**: [Snakefile](file://e:\github\Bioit-CEITEC\demultiplexing_short_read\Snakefile)
- Creates sample table from config using `get_panda_sample_tab_from_config()`
- Expands output FASTQ files based on read count (R1, R2, R3/UMI)
- Includes platform-specific demultiplexing rules

**Rules** ([rules/](file://e:\github\Bioit-CEITEC\demultiplexing_short_read\rules))
- [demultiplexing.smk](file://e:\github\Bioit-CEITEC\demultiplexing_short_read\rules\demultiplexing.smk): Platform-specific demultiplexing (Aviti `Bases2Fastq`, MGI `calDemux`, Illumina `bcl2fastq`) + `fastq_mv` for reorganizing output
- [merge.smk]: Combines FASTQ files from multiple libraries into single output
- [fastqc.smk]: Quality control with FastQC and MultiQC

**Wrappers** ([wrappers/](file://e:\github\Bioit-CEITEC\demultiplexing_short_read\wrappers))
Each wrapper contains a `script.py` and `env.yaml` for conda environment:
- `aviti_bases2fastq/script.py`: Aviti Bases to FASTQ conversion
- `mgi_calDemux/script.py`: MGI demultiplexing via SplitBarcode
- `bcl2fastq/script.py`: Illumina bcl2fastq execution
- `fastq_mv/script.py`: Moves/cuts FASTQ files to final locations
- `stats_copy/script.py`: Extracts demultiplexing statistics

## Output Format

Final FASTQ files: `{library}/raw_fastq/{sample_name}_R{read_num}.fastq.gz`
- R1: Forward reads
- R2: Reverse reads (if paired-end)
- R3/UMI: UMI reads (Aviti platform with UmiMask configured)

## Running the Workflow

```bash
snakemake --snakefile Snakefile --configfile config.json --cores <n>
```

For merged output mode, set `"merged": true` in config.json.

## Platform-Specific Details

**Aviti**: Uses `Bases2Fastq` tool with base masks (R1FastQMask, R2FastQMask, UmiMask) defined per library
**MGI**: Uses `SplitBarcode` tool with lane-based processing
**Illumina**: Uses `bcl2fastq` with generated sample sheets

## Common Patterns

- Sample table generation uses pandas DataFrame with library/sample metadata
- Wildcard constraints are dynamically generated from sample_tab
- Demux settings group samples with identical demultiplexing parameters
- Lane splitting creates separate processing paths per lane (L001, L002, etc.)
