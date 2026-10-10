# Multi-omics Sequencing Data Processing Pipelines

Core command-line workflows for paired-end whole-genome sequencing (WGS), bisulfite sequencing (5mC workflow), 5hmC enrichment sequencing, and RNA-seq. This README consolidates four supplied workflows into reusable examples with English comments and generic filenames.

These are command templates, not software defaults or a fully automated workflow. Replace the input, reference, annotation, and sample identifiers before running. Commands have been reviewed for consistency but have not been executed on sequencing data.

## Pipeline overview

| Workflow | Quality control / trimming | Alignment | Duplicate removal | Main output |
| --- | --- | --- | --- | --- |
| WGS | fastp | BWA-MEM2 | GATK MarkDuplicates | Coordinate-sorted, deduplicated genomic BAM |
| Bisulfite sequencing / 5mC | FastQC, MultiQC, Trim Galore | Bismark | deduplicate_bismark | Deduplicated BAM, cytosine methylation tables, regional summaries |
| 5hmC enrichment sequencing | FastQC, MultiQC, Trimmomatic | Bowtie2 | Picard MarkDuplicates | Filtered, coordinate-sorted, deduplicated genomic BAM |
| RNA-seq | FastQC, MultiQC, Trimmomatic | STAR | Picard MarkDuplicates | Genomic BAM, transcriptome BAM, gene and transcript-grouped count tables |

The third workflow is labeled as 5hmC enrichment sequencing based on the supplied project context. Its commands perform standard DNA read alignment and preprocessing; assay-specific enrichment quantification is outside this README.

## Requirements and conventions

- Use Bash, with the relevant executables available in `PATH`.
- Install the tools listed in the overview. Bismark also requires the alignment backend supported by the installed version; legacy Trim Galore installations require Cutadapt. Trimmomatic requires Java.
- Supply matching paired-end FASTQ files: `input_R1.fq.gz` and `input_R2.fq.gz`.
- Supply an uncompressed reference FASTA named `reference.fa`. Use the same assembly and chromosome naming throughout alignment, methylation analysis, and annotation.
- Supply `annotation.gtf` for STAR and featureCounts, and `TruSeq3-PE.fa` for the Trimmomatic examples. Replace `trimmomatic.jar` with the installed JAR filename or path.
- Run each sample and each workflow in a separate working directory to avoid overwriting generic output names.
- Build each tool's reference index once per reference and relevant configuration. Bismark and STAR use index directories; these are required tool inputs, not project-specific paths.
- Set `RGID`, `RGLB`, `RGPU`, and `RGSM` to the actual read group, library, platform unit, and sample identifiers. Use distinct read groups for distinct sequencing units when appropriate.

Before running a block, enable pipeline error detection:

```bash
set -euo pipefail
```

Record software versions and the reference/annotation release for reproducibility. Thread counts largely follow the supplied workflows and should be adjusted to the available hardware. `samtools sort -m 10G` sets an approximate memory limit per thread, so the examples can require substantial RAM. A pipeline containing multiple threaded tools uses resources concurrently.

## 1. Whole-genome sequencing

### Quality filtering, alignment, sorting, and duplicate removal

```bash
# Filter and trim paired-end reads
fastp -w 16 \
    --qualified_quality_phred 30 \
    --length_required 36 \
    --cut_mean_quality 30 --cut_front --cut_tail \
    --trim_front1 5 --trim_tail1 3 \
    --trim_front2 5 --trim_tail2 3 \
    -i input_R1.fq.gz -I input_R2.fq.gz \
    -o clean_R1.fq.gz -O clean_R2.fq.gz \
    -h fastp.html -j fastp.json

# Aggregate quality-control reports
multiqc .

# Build the reference index once
bwa-mem2 index reference.fa

# Align paired-end reads
bwa-mem2 mem -t 16 reference.fa clean_R1.fq.gz clean_R2.fq.gz \
    > aligned.sam

# Sort genomic alignments by coordinates
samtools sort -@ 20 -l 6 -m 10G aligned.sam -o sorted.bam

# Add read group information
picard AddOrReplaceReadGroups \
    I=sorted.bam O=sorted_RG.bam \
    RGID=sample RGLB=library RGPL=ILLUMINA RGPU=unit RGSM=sample

# Remove duplicate reads and write duplication metrics
gatk MarkDuplicates \
    I=sorted_RG.bam O=deduplicated.bam \
    M=duplicate_metrics.txt \
    REMOVE_DUPLICATES=true
```

The main outputs are cleaned FASTQ files, fastp reports, `deduplicated.bam`, and `duplicate_metrics.txt`. Variant calling and base-quality recalibration are outside this workflow.

## 2. Bisulfite sequencing / 5mC workflow

### Quality control, trimming, alignment, and duplicate removal

Place `reference.fa` in `bismark_genome` before preparing the index. This example assumes a directional paired-end bisulfite library; adapt alignment settings if the library preparation differs.

```bash
# Assess raw read quality
fastqc -t 60 input_R1.fq.gz input_R2.fq.gz
multiqc .

# Trim adapters and low-quality bases
# Eight trimming cores are used as a portable example
trim_galore --cores 8 --illumina --paired --length 36 --gzip \
    --clip_R1 5 --clip_R2 5 \
    --three_prime_clip_R1 3 --three_prime_clip_R2 3 \
    --quality 30 \
    input_R1.fq.gz input_R2.fq.gz

# Prepare the bisulfite reference index once
bismark_genome_preparation --parallel 16 bismark_genome

# Align trimmed paired-end reads
bismark -p 2 --gzip --score_min L,0,-0.2 \
    --genome bismark_genome \
    -1 input_R1_val_1.fq.gz \
    -2 input_R2_val_2.fq.gz

# Remove duplicates before coordinate sorting
deduplicate_bismark --paired --bam \
    input_R1_val_1_bismark_bt2_pe.bam

# Sort the deduplicated BAM for downstream analysis
samtools sort -@ 32 -l 6 -m 10G \
    input_R1_val_1_bismark_bt2_pe.deduplicated.bam \
    -o sorted.bam

# Generate alignment statistics and aggregate reports
samtools stats sorted.bam > alignment_stats.txt
multiqc .

# Assess the processed BAM
fastqc -t 60 sorted.bam
multiqc .
```

Trim Galore normally produces paired outputs ending in `_val_1.fq.gz` and `_val_2.fq.gz`. The BAM names above follow the classic Bismark naming convention; use the filenames generated by the installed release. Deduplicate the original paired Bismark output before coordinate sorting so mate records remain suitable for its duplicate-removal procedure.

### Methylation calculation and feature summaries

The following commands retain the BatMeth2 interface from the supplied workflow. BatMeth2 releases differ in indexing and annotation subcommands: some document `build_index` and `annoation` instead of `index` and `methyGff`. Confirm the interface and compatibility of externally aligned Bismark BAM files with the installed release before running this section.

```bash
# Prepare the BatMeth2 reference data once
BatMeth2 index -g reference.fa

# Calculate cytosine methylation from coordinate-sorted alignments
BatMeth2 calmeth \
    -Q 30 --coverage 4 -nC 1 \
    --Regions 200 --step 50000 \
    --genome reference.fa \
    --binput sorted.bam \
    --methratio sample

# Summarize gene bodies and promoters with 2 kb flanking regions
BatMeth2 methyGff \
    --coverage 4 -nC 1 \
    --genome reference.fa \
    -b5 genes.bed5 \
    --body --promoter --GENE \
    --distance 2000 --step 0.01 \
    --methratio sample.methratio.txt \
    --out sample.genes

# Summarize transcript bodies and promoters with 1 kb flanking regions
BatMeth2 methyGff \
    --coverage 4 -nC 1 \
    --genome reference.fa \
    -b5 transcripts.bed5 \
    --body --promoter --GENE \
    --distance 1000 --step 0.01 \
    --methratio sample.methratio.txt \
    --out sample.transcripts
```

For BatMeth2 `-b5`, supply five tab-separated columns in this order: chromosome, start, end, feature ID, strand. Use the coordinate convention required by the installed release; when converting a 1-based inclusive GTF to standard BED, subtract one from the start and leave the end unchanged. Parse GTF attributes by name rather than assuming a fixed whitespace column contains the gene ID.

The main outputs include `sample.methratio.txt` and gene/transcript regional summaries. The assay chemistry determines how the methylation signal should be interpreted; this preprocessing workflow alone does not establish modification specificity.

## 3. 5hmC enrichment sequencing

### Quality control, trimming, alignment filtering, and duplicate removal

```bash
# Assess raw read quality
fastqc -t 60 input_R1.fq.gz input_R2.fq.gz
multiqc .

# Trim adapters and filter low-quality reads
java -jar trimmomatic.jar PE -threads 40 \
    input_R1.fq.gz input_R2.fq.gz \
    clean_R1_paired.fastq.gz clean_R1_unpaired.fastq.gz \
    clean_R2_paired.fastq.gz clean_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:5 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:36

# Build the reference index once
bowtie2-build --threads 60 reference.fa reference_index

# Align paired reads, retain mapped alignments with MAPQ >= 10, and sort
bowtie2 --threads 20 --end-to-end --very-sensitive --no-unal \
    -x reference_index \
    -1 clean_R1_paired.fastq.gz \
    -2 clean_R2_paired.fastq.gz \
    | samtools view -b -q 10 -F 4 \
    | samtools sort -@ 20 -l 5 -m 10G -o sorted.bam

# Add read group information
picard AddOrReplaceReadGroups \
    I=sorted.bam O=sorted_RG.bam \
    RGID=sample RGLB=library RGPL=ILLUMINA RGPU=unit RGSM=sample

# Remove duplicate reads and write duplication metrics
picard MarkDuplicates \
    I=sorted_RG.bam O=deduplicated.bam \
    M=duplicate_metrics.txt \
    REMOVE_DUPLICATES=true
```

Only paired Trimmomatic outputs are aligned here. The `samtools view` filter operates on individual alignment records and does not impose a proper-pair requirement. Coordinate sorting is already performed in the alignment pipeline, so a second sorting pass is unnecessary.

The main output is `deduplicated.bam`. Peak calling, enrichment normalization, and differential 5hmC analysis are outside this workflow.

## 4. RNA-seq

### Quality control, trimming, and STAR alignment

Create `star_genome` before building the index. STAR genome generation and alignment can require substantial RAM for a human reference.

```bash
# Assess raw read quality
fastqc -t 60 input_R1.fq.gz input_R2.fq.gz
multiqc .

# Trim adapters and filter low-quality reads
java -jar trimmomatic.jar PE -threads 40 \
    input_R1.fq.gz input_R2.fq.gz \
    clean_R1_paired.fastq.gz clean_R1_unpaired.fastq.gz \
    clean_R2_paired.fastq.gz clean_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:5 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:36

# Build the genome index with matching gene annotation once
STAR --runMode genomeGenerate \
    --runThreadN 40 \
    --genomeDir star_genome \
    --genomeFastaFiles reference.fa \
    --sjdbGTFfile annotation.gtf

# Generate genomic alignments, transcriptome alignments, and STAR gene counts
STAR --runMode alignReads \
    --runThreadN 40 \
    --readFilesCommand zcat \
    --genomeDir star_genome \
    --readFilesIn clean_R1_paired.fastq.gz clean_R2_paired.fastq.gz \
    --quantMode TranscriptomeSAM GeneCounts \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix sample_

# Sort genomic alignments for annotation-based counting
samtools sort -@ 20 -l 6 -m 10G \
    sample_Aligned.out.bam -o sorted.bam

# Optionally sort transcriptome alignments for a compatible downstream tool
samtools sort -@ 20 -l 6 -m 10G \
    sample_Aligned.toTranscriptome.out.bam -o transcriptome_sorted.bam
```

Use the genomic BAM with a genomic GTF in featureCounts. The transcriptome BAM has transcript coordinates and should only be used with a downstream tool that explicitly accepts that format and ordering.

### Read groups and duplicate removal

This section reproduces the supplied duplicate-removal workflow. Keep the pre-removal BAM so the counting strategy can be chosen for the experiment.

```bash
# Add read group information to genomic alignments
picard AddOrReplaceReadGroups \
    I=sorted.bam O=sorted_RG.bam \
    RGID=sample RGLB=library RGPL=ILLUMINA RGPU=unit RGSM=sample

# Remove duplicate reads and write duplication metrics
picard MarkDuplicates \
    I=sorted_RG.bam O=deduplicated.bam \
    M=duplicate_metrics.txt \
    REMOVE_DUPLICATES=true
```

### Gene and transcript-grouped counts

The examples below count individual reads, preserving the supplied `-p` commands for Subread 2.0.2 and later. In those releases, `-p` declares paired-end input; add `--countReadPairs` to count fragments instead. Earlier releases gave `-p` different counting behavior, so record the version and counting unit.

```bash
# Count exon-overlapping reads grouped by gene ID
# Add --countReadPairs if fragment counts are intended
featureCounts -T 32 -t exon -g gene_id -p \
    -a annotation.gtf -F GTF \
    -o gene_counts.txt \
    deduplicated.bam

# Count exon-overlapping reads grouped by transcript ID
# Add --countReadPairs if fragment counts are intended
featureCounts -T 36 -t exon -g transcript_id -p \
    -a annotation.gtf -F GTF \
    -o transcript_counts.txt \
    deduplicated.bam

# Optionally create a copy without metadata comments
sed '/^#/d' transcript_counts.txt > transcript_counts_no_comments.txt
```

Choose the strandedness setting (`-s 0`, `-s 1`, or `-s 2`) according to the library preparation; the examples use the default unstranded setting. Transcript-grouped featureCounts output is an overlap-based count table rather than an isoform abundance estimate, and shared exons can cause ambiguous assignments. STAR's `ReadsPerGene.out.tab` is produced before the later duplicate-removal step and should not be assumed to match these tables.

The main outputs are genomic and transcriptome BAM files, `gene_counts.txt`, `transcript_counts.txt`, and their assignment summaries. Differential expression analysis is outside this workflow.

## Parameter interpretation

| Tool / option | Meaning in these examples |
| --- | --- |
| fastp `--qualified_quality_phred 30` | Threshold used to classify qualified bases; it does not require every retained base to have Q30 |
| fastp `--cut_mean_quality 30` | Mean quality threshold for enabled end-cutting operations |
| fastp / Trim Galore end clipping | Remove 5 bases from the 5-prime end and 3 bases from the 3-prime end of each mate |
| Trimmomatic `LEADING:5 TRAILING:3` | Quality thresholds for trimming end bases, not fixed numbers of bases to remove |
| Trimmomatic `SLIDINGWINDOW:5:15` | Trim when the mean quality in a five-base window falls below 15 |
| Minimum length 36 | Minimum read length used by the trimming examples |
| `samtools view -q 10 -F 4` | Retain mapped alignment records with MAPQ at least 10 |
| `REMOVE_DUPLICATES=true` | Remove duplicate records rather than retaining marked records |
| BatMeth2 `--coverage 4 -nC 1` | Coverage and cytosine-count thresholds used in the supplied methylation workflow |

## Reuse and validation

Run stages sequentially and inspect each tool's report before proceeding. Avoid launching an unlimited number of background jobs. Review trimming retention, alignment rate, duplication metrics, and the featureCounts assignment summaries where applicable.

All four examples retain the supplied duplicate-removal approach. Decide whether removal is appropriate for the assay; coordinate-based duplicate removal does not substitute for UMI-aware processing. Preserve intermediate BAM files until the analysis choices and outputs have been checked.

Optionally index final coordinate-sorted genomic BAM files for downstream random access:

```bash
# Use deduplicated.bam for WGS, 5hmC, or the RNA-seq example
samtools index deduplicated.bam

# Use sorted.bam for the bisulfite example
samtools index sorted.bam
```

The original directory creation, symlinks, batch-script generation, background jobs, and sample-specific renaming have been omitted. The duplicate FastQC thread option and redundant Bowtie2 sorting pass have been consolidated. The optional count-table cleanup writes a new file instead of editing the original output.

## Tool documentation

- [fastp](https://github.com/OpenGene/fastp)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](https://docs.seqera.io/multiqc)
- [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2)
- [SAMtools](https://www.htslib.org/doc/samtools.html)
- [Picard](https://broadinstitute.github.io/picard/)
- [GATK MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard)
- [Trim Galore](https://github.com/FelixKrueger/TrimGalore)
- [Bismark duplicate removal](https://felixkrueger.github.io/Bismark/usage/deduplication/)
- [BatMeth2](https://github.com/GuoliangLi-HZAU/BatMeth2)
- [BatMeth2 methylation calculation](https://batmeth2-docs.readthedocs.io/en/latest/function/Calmeth.html)
- [Trimmomatic](https://github.com/usadellab/Trimmomatic)
- [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [STAR](https://github.com/alexdobin/STAR)
- [featureCounts](https://subread.sourceforge.net/featureCounts.html)

Consult the documentation for the exact installed release. No complete environment specification or software-version manifest was supplied with the original workflows.
