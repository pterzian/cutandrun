# CUT&RUN Data Processing Pipeline - using Unique Molecular Identifiers (UMI)

This repository provides a pipeline for processing CUT&RUN data, encompassing steps for trimming reads, aligning them to both primary and spike-in genomes, deduplicating reads, and calling peaks.

> [!WARNING]
> You are on the "umi-tools" branch of this repo which is still under development.


## Pipeline Overview

1. **Trim Reads:** Remove adapter sequences and low-quality bases from raw sequencing reads.
2. **Align Reads to Primary Genome:** Map trimmed reads to the primary reference genome.
3. **Align Reads to Spike-in Genome:** Map trimmed reads to the spike-in reference genome for normalization purposes.
4. **Deduplicate Reads using UMI-tools :** Eliminate PCR duplicates to ensure data accuracy.
 

## Getting Started

Ensure that the following tools are installed and accessible in your system's PATH:

- [Fastp](https://github.com/OpenGene/fastp): For trimming adapter sequences
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml): For aligning reads to reference genomes
- [Samtools](http://www.htslib.org/): For handling SAM/BAM files
- [UMI-tools](https://github.com/CGATOxford/UMI-tools): For deduplication of reads
- [Deeptools](https://github.com/deeptools/deepTools/): For downstream analysis <mark>(not needed to run the pipeline)</mark>

Or, you can use the provided conda yaml recipe to create a suited environment (it will install both UMI-tools and Picard tools).

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

Sample sheets can be built using `./build_sample_sheet.sh /path/to/raw_reads_folders`

> [!NOTE] 
> When sequencing required multiple runs to reach the expected sequencing depth, each run much have its own line in the sample_sheet. **Reads with the same sample name** will be merged at the start of the pipeline.*


### Expected output of the pipeline

* One folder per sample with the following structure :

```
|-WT1_H3k27ac
 |-mapped_spikein
 | |-WT1_H3k27ac_spikein_flagstat.txt
 | |-WT1_H3k27ac_spikein_aligned.bam
 |-deduplicated_umi
 | |-primary
 | | |-WT1_H3k27ac_primary_dedup_UMI.bam
 | | |-WT1_H3k27ac_primary_dedup_UMI_dedup.log
 | |-spikein
 | | |-WT1_H3k27ac_spikein_dedup_UMI.bam
 | | |-WT1_H3k27ac_spikein_dedup_UMI_dedup.log
 |-merged_fastq
 | | |-WT1_H3k27ac_merged_R1.fq.gz
 | | |-WT1_H3k27ac_merged_R2.fq.gz
 |-trimmed
 | |-WT1_H3k27ac_trimmed_R2.fq.gz
 | |-fastp_report.json
 | |-WT1_H3k27ac_trimmed_R1.fq.gz
 | |-fastp_report.html
 |-mapped_primary
 | |-WT1_H3k27ac_primary_flagstat.txt
 | |-WT1_H3k27ac_primary_aligned.bam
 |-umi
 | |-WT1_H3k27ac_umi_R2.fq.gz
 | |-WT1_H3k27ac_umi_extract.log
 | |-WT1_H3k27ac_umi_R1.fq.gz.fq.gz
```

* A file called `flagstat_summary.tsv` collecting metrics about the number of reads at the different steps of the pipeline that can be usefull for downstram spike in normalization.

* A `timing_report.tsv` file with the time spend for each steps of the pipeline


### Downstream analysis

**Correlation between samples**

The script `run_plot_cor.sh` was designed to generate scatterplots, heatmaps, and PCA plots using `deeptools` library. 

```
Usage: run_plot_cor.sh <study_name> <bam_dir> <bin_size> [bam_pattern]
```

Here is an example :

```
run_plot_cor.sh test_run ../all_WT_bam/ 200  "S1_WT2_D_0.001*.bam"
```

> [!WARNING]
> The first step of `run_plot_cor.sh` is the generation of a `multiBamSummary` matrix that needs all bam to be indexed. Because the last step of the pipeline is the deduplication, output bams are not sorted or indexed. This supplementary step can be done using `samtools sort` and `samtools index` or `sambamba sort` that combine both sorting and indexing.


**Bam to bigwig file**

The script `generate_bamcoverage_command.sh` was written to generates list of bamcoverage commands that can be launch to generate bigwig files.

```
Usage: ./generate_bamcoverage_command.sh <bam_dir> <bam_pattern> <normalize_method: CPM|RPKM> <blacklist.bed> <stranded: stranded|unstranded> <output_script.sh> [output_dir] [--binSize N] [--nproc N]

Parameters:
  bam_dir          Directory containing BAM files.
  bam_pattern      Pattern to match BAM files (e.g., '*.bam').
  normalize_method Normalization method: CPM or RPKM.
  blacklist.bed    Path to blacklist BED file.
  stranded         'stranded' to generate stranded bigWigs (forward & reverse), 'unstranded' for unstranded.
  output_script.sh Output bash script file name to write commands.
  output_dir       (Optional) Directory to save bigWig files. Defaults to current directory.

Optional flags:
  --binSize N      Bin size for bamCoverage (default: 5).
  --nproc N        Number of processors/threads for bamCoverage (default: 2).

Example:
  ./generate_bamcoverage_command.sh /data/bams '*.bam' CPM /refs/blacklist.bed stranded bigwig_commands.sh /results --binSize 10 --nproc 4
```

> [!IMPORTANT]  
> Make sure to quote the pattern to match BAM files.

The output of this script is file containing all `bamcoverage` command lines.