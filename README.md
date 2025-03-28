# CUT&RUN Data Processing Pipeline

This repository provides a pipeline for processing CUT&RUN data, encompassing steps for trimming reads, aligning them to both primary and spike-in genomes, deduplicating reads, and calling peaks.

## Pipeline Overview

1. **Trim Reads:** Remove adapter sequences and low-quality bases from raw sequencing reads.
2. **Align Reads to Primary Genome:** Map trimmed reads to the primary reference genome.
3. **Align Reads to Spike-in Genome:** Map trimmed reads to the spike-in reference genome for normalization purposes.
4. **Deduplicate Reads:** Eliminate PCR duplicates to ensure data accuracy.
5. **Call Peaks:** Identify regions of significant enrichment, indicating protein-DNA interactions.

## Getting Started

### Prerequisites

Ensure that the following tools are installed and accessible in your system's PATH:

- [Fastp](https://github.com/OpenGene/fastp): For trimming adapter sequences.
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml): For aligning reads to reference genomes.
- [Samtools](http://www.htslib.org/): For handling SAM/BAM files.
- [Picard](http://broadinstitute.github.io/picard/): For deduplication of reads.
- [MACS2](https://github.com/macs3-project/MACS): For peak calling.

## Use the pipeline

The pipeline can be launched from the following command-line :

```
./run_cutandrun_2.sh samplesheet.csv path/to/bowtie2-index-target path/to/bowtie2-index-spikein
```

The `samplesheet.csv` is comma-separated table with 4 fields such as : 

```
SampleName,AssayName,FASTQ_R1,FASTQ_R2
WT1_IgG,cutandrun,path/to/WT1_IgG_1.fq.gz,path/to/WT1_IgG_2.fq.gz
WT1_H3k27ac,cutandrun,/path/to/WT1_H3k27ac_xxx_1.fq.gz,/path/to/WT1_H3k27ac_xxx_2.fq.gz
WT1_H3k27ac,cutandrun,/path/to/WT1_H3k27ac_xyx_1.fq.gz,/path/to/WT1_H3k27ac_xyx_2.fq.gz
```

### Easy building of the sample sheet

Sample sheets can be build using `./build_sample_sheet.sh /path/to/raw_reads_folders`

*Notes : When sequencing required multiple run to reach the expected sequencing depth, each run much have its own line in the sample_sheet. Reads with the same sample name will be merged at the start of the pipeline.*

 
