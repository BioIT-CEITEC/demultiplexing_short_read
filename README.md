# Demultiplexing Short Read Workflow

A Snakemake-based workflow for demultiplexing short-read sequencing data from multiple platforms (Illumina, Aviti, MGI). This workflow converts raw sequencing output into sample-specific FASTQ files ready for downstream analysis.

## Features

- **Multi-platform support**: Illumina (bcl2fastq), Aviti (Bases2Fastq), and MGI (SplitBarcode)
- **Flexible sample configuration**: Support for multiple libraries and samples per run
- **Lane splitting**: Optional processing by sequencing lane
- **UMI support**: Handle UMI reads (Aviti platform)
- **Quality control**: Integrated FastQC and MultiQC reporting
- **Read length trimming**: Optionally cut reads to library-specific lengths

## Workflow Overview

```
Raw Sequencing Data → Demultiplexing → FASTQ Organization → Quality Control
     (BCL/Bases)         (Platform)        (raw_fastq/)        (FastQC)
```

### Pipeline Steps

1. **Sample Sheet Generation**: Creates platform-specific sample sheets from configuration
2. **Demultiplexing**: Platform-specific conversion of raw data to FASTQ
3. **FASTQ Reorganization**: Moves and optionally trims FASTQ files to final locations
4. **Statistics Collection**: Extracts demultiplexing metrics (read counts, etc.)
5. **Quality Control**: Runs FastQC on all FASTQ files and generates MultiQC reports

## Configuration

Create a `config.json` file with the following structure:

```json
{
  "run_dir": "/path/to/sequencing/run",
  "run_date": "2024-01-15",
  "run_name": "MyRun",
  "run_sequencer_type": "AVITI",
  "run_forward_read_length": 52,
  "run_reverse_read_length": 52,
  "run_lane_splitting": null,
  "globalTmpdPath": "/tmp",
  "globalResources": "/path/to/resources",
  "libraries": {
    "Library1": {
      "samples": {
        "Sample1": {
          "sample_name": "Sample1",
          "i7_sequence": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
          "i5_sequence": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
        }
      },
      "lib_forward_read_length": 52,
      "lib_reverse_read_length": 52,
      "cut_reads_to_lib_size": false,
      "barcode_mismatches": 0,
      "AVITI_R1FastQMask": "YYYYYYYYYYYYYYYYYYYYYYYYYYYYYY",
      "AVITI_R2FastQMask": "YYYYYYYYYYYYYYYYYYYYYYYYYYYYYY",
      "AVITI_UmiMask": ""
    }
  }
}
```

### Sequencer Types

| Type | Description | Demultiplexing Tool |
|------|-------------|---------------------|
| `ILLUMINA` | Illumina sequencers | bcl2fastq |
| `AVITI` | Element Biosciences AVITI | Bases2Fastq |
| `MGI` | MGI sequencers | SplitBarcode |

### Library Configuration Options

**Common options:**
- `samples`: Dictionary of samples with barcode sequences
- `lib_forward_read_length`: Target forward read length
- `lib_reverse_read_length`: Target reverse read length
- `cut_reads_to_lib_size`: Whether to trim reads to library size
- `barcode_mismatches`: Allowed barcode mismatches

**Platform-specific options:**
- **Aviti**: `AVITI_R1FastQMask`, `AVITI_R2FastQMask`, `AVITI_I1Mask`, `AVITI_I2Mask`, `AVITI_UmiMask`
- **MGI**: `MGI_filter_param`

## Usage

### Basic Execution

```bash
snakemake --cores <num_cores> --configfile config.json
```

### With Dry Run (validate workflow)

```bash
snakemake --dry-run --cores 1 --configfile config.json
```

### Generate DAG Visualization

```bash
snakemake --dag | dot -Tpdf > dag.pdf
```

### Merge Mode

To merge FASTQ files from multiple libraries into a single output directory, set `"merged": true` in your config.json and specify `library_output`.

## Output Structure

```
<library>/
├── raw_fastq/
│   ├── <sample_name>_R1.fastq.gz    # Forward reads
│   ├── <sample_name>_R2.fastq.gz    # Reverse reads (if paired-end)
│   └── <sample_name>_UMI.fastq.gz   # UMI reads (if configured)
├── qc_reports/
│   ├── <sample_name>/
│   │   └── raw_fastqc/
│   │       └── <sample_name>_R*_fastqc.html
│   └── raw_fastq_multiqc.html
└── sequencing_run_info/
    └── samplesNumberReads.json
```

## Dependencies

The workflow uses conda environments for each wrapper. Required software:
- **Snakemake** >= 5.18.0
- **Python** 3.10+
- **Pandas** (for sample table generation)
- **Cutadapt** (for read trimming)
- Platform-specific tools (configured via conda environments)

## File Organization

```
.
├── Snakefile                    # Main workflow definition
├── config.json                  # Runtime configuration (not in repo)
├── AGENT.md                     # Claude Code instructions
├── README.md                    # This file
├── rules/
│   ├── demultiplexing.smk       # Demultiplexing rules per platform
│   ├── fastqc.smk               # Quality control rules
│   └── merge.smk                # Merge mode rules
└── wrappers/
    ├── aviti_bases2fastq/       # Aviti Bases2Fastq wrapper
    ├── bcl2fastq/               # Illumina bcl2fastq wrapper
    ├── mgi_calDemux/            # MGI calDemux wrapper
    ├── fastq_mv/                # FASTQ reorganization wrapper
    ├── stats_copy/              # Statistics extraction wrapper
    └── ...                      # Other wrappers
```

## Troubleshooting

### Common Issues

1. **"No rule to make target"**: Ensure your `config.json` is properly configured and all required paths exist.

2. **Missing conda environments**: Run `snakemake --use-conda --cores 1` to create missing environments.

3. **Demultiplexing errors**: Check that barcode sequences in config match the actual library prep and that base masks are correctly specified for Aviti runs.

## License

This project is part of the CEITEC Bioinformatics workflows collection.
