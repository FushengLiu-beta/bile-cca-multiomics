# Multi-omics Sequencing Data Processing Pipelines

Core command-line workflows for paired-end whole-genome sequencing (WGS), bisulfite sequencing (5mC workflow), 5hmC enrichment sequencing, and RNA-seq. This README consolidates four supplied workflows into reusable examples with English comments and generic filenames.

These are command templates, not software defaults or a fully automated workflow. Replace the input, reference, annotation, and sample identifiers before running. Commands have been reviewed for consistency.

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



## 5. FOREST: Frequency-based tumor classification and out-of-fold predictions

This section adds a standalone R workflow after the sequencing workflows above. It accepts a sample-by-feature CSV, repeats stratified training/validation splits, selects features by their training-only inclusion frequency, and writes validation predictions to the current working directory.

This version runs directly with `train_fraction = 0.67`, `TOP_n = 10L`, `alpha_grid = seq(0.5, 1, by = 0.05)`, and `lambda_ratio_min = 0.001`. No preliminary parameter-screening run or prior analysis object is required. In each outer split, feature frequencies are pooled across alpha/lambda paths fitted on the complete training set. Training CV then selects one final alpha/lambda pair, and the outer validation set evaluates that model.

**FOREST** stands for **F**eature-frequency **O**ut-of-fold **R**esampling with **E**lastic-net **S**election and **T**raining. The name describes the workflow's feature-frequency selection, repeated holdout evaluation, and elastic-net training. Its script filename is `FOREST.R`.

The supplied files are `FOREST.R`, `input.csv`, and this extended `README.md`. The R script includes all modeling and export functions; it does not need `frequency_p.R` or another project script.

### Input format

Use a UTF-8, comma-separated file with a header and one row per independent sample.

| Position | Content | Requirements |
| --- | --- | --- |
| Column 1 | Sample ID | Unique, nonempty text; leading zeroes are preserved |
| Column 2 | Sample group | Exactly `Tumor` or `Control`; both groups must be present |
| Columns 3 onward | Feature expression | At least two numeric feature columns with unique, nonempty names |

The first two headers are internally renamed to `ID` and `Group`; feature names are preserved. Do not include an exported row-index column. Metadata such as age, batch, and sex should not be included among the expression columns unless intentionally treated as model features.

Missing expression values may be empty, `NA`, or `NaN`. The script rejects invalid numeric text, infinite values, empty/duplicate IDs, missing groups, and duplicate headers. It trims surrounding whitespace from IDs and groups. Feature columns cannot be named `ID` or `Group`.

Provide feature values that have already undergone the assay-appropriate processing needed to obtain an expression matrix. This classification script does not process FASTQ/BAM files or perform RNA-seq count normalization. It estimates its imputation and standardization parameters using each training set.

The example `input.csv` contains **80 simulated samples**, 40 in each group, and 30 features. Twenty features have engineered group shifts, nine are noise features, and one is constant; two expression cells are missing. Values use Gaussian noise generated with a fixed Mulberry32 seed of 20261010 and are rounded to six decimal places. This is a format and workflow example; its AUC does not represent biological performance.

### Run in an R session

An R environment is required. Put the script and your `input.csv` in the current working directory, then run:

```r
getwd()
n_repeats <- 100L
frequency_config <- list()        # Use the script's fixed defaults
source("FOREST.R", encoding = "UTF-8")
```

The script runs immediately. It checks for `glmnet` and `pROC` and, by default, attempts to install missing packages from CRAN. Package installation needs internet access and a writable R library; Linux source installations may also need the relevant compiler toolchain. To install dependencies explicitly once:

```r
install.packages(c("glmnet", "pROC"), repos = "https://cloud.r-project.org")
```

Results are written to `getwd()` by default. Source execution also creates `frequency_run`, `results_df`, `model_df`, `auc_df`, `curve_df`, `oof_predictions`, `oof_predictions_summary`, and `indexfilterresult` in the sourcing environment.

### Configure parameters on the current page

Set `n_repeats <- 100L` on the current page before sourcing the script. It controls the number of outer training/validation splits. Its default is 100. Feature-frequency estimation uses the complete outer training set for each alpha/lambda path.

For other settings, define `frequency_config` before sourcing the script. Only supplied entries override defaults. An explicit `frequency_config$n_repeats` takes precedence over the page-level `n_repeats`; otherwise the page-level value is used. Each formal run uses one fixed training fraction, TOP_n, lambda lower ratio, and an explicitly specified alpha grid.

The defaults below are the starting settings for this standalone workflow. Parameters are specified directly on the current page; they are not estimated by a separate parameter-screening stage.

```r
n_repeats <- 100L
frequency_config <- list(
  input_file = "input.csv",
  output_dir = ".",                # Current working directory
  train_fraction = 0.67,
  TOP_n = 10L,
  alpha_grid = seq(0.5, 1, by = 0.05), # 11 alpha values
  lambda_ratio_min = 0.001,
  n_repeats = n_repeats,            # Outer training/validation splits
  n_lambda = 50L,
  pvalue_cutoff = 0.05,
  seed = 123L,
  overwrite = FALSE
)
source("FOREST.R", encoding = "UTF-8")
```

To disable differential p-value filtering, include `pvalue_cutoff = NULL` in the configuration **before** the run. NULL is retained as an intentional override. All features then enter the candidate pool, subject only to training-set numerical usability checks.

| Setting | Default | Meaning |
| --- | --- | --- |
| `positive_class` | `"Tumor"` | Positive class for this CSV entry point; probabilities always refer to Tumor |
| `train_fraction` | 0.67 | Per-class training count is `floor(class_count * train_fraction)`. A value of 0.67 splits each class approximately 2:1 between training and validation. |
| `TOP_n` | 10 | Number of highest positive-frequency features supplied to each final model |
| `alpha_grid` | `seq(0.5, 1, by = 0.05)` | 11 elastic-net mixing values used for pooled feature frequencies and final training CV |
| `lambda_ratio_min` | 0.001 | Lower end of a log-spaced lambda/lambda_max path; upper ratio is 1 |
| `n_repeats` | 100 | Number of outer splits and final model attempts |
| `n_lambda` | 50 | Requested points in each lambda path |
| `pvalue_cutoff` | 0.05 | Training-only raw p-value threshold; NULL disables this step |
| `diff_method` | `"wilcox"` | Two-sided Wilcoxon test; `"welch"` uses Welch's t-test |
| `cv_nfolds` | 5 | Final training-only CV folds; reduced when a training class is smaller |
| `plot_iter` | 1 | Split used by the optional ROC plotting function; all repeat details are saved |
| `fpr_grid` | 101 points from 0 to 1 | Grid used for interpolated ROC curve export |
| `seed` | 123 | Controls outer splits, coverage repair, and training CV fold assignment |
| `coef_tol` | 1e-8 | Absolute coefficient threshold for counting feature inclusion |
| `fit_maxit` | 1000000 | Maximum glmnet iterations |
| `ensure_oof_coverage` | TRUE | Repair the planned splits so every sample appears in validation |
| `require_complete_oof` | TRUE | Raise an error after saving diagnostics if any sample has no valid OOF prediction |
| `install_missing_packages` | TRUE | Attempt to install missing glmnet/pROC packages |
| `save_models` | TRUE | Save the complete result and models as `frequency_analysis.rds` |
| `overwrite` | FALSE | Refuse to overwrite existing named outputs |

The default alpha grid is `0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00`. Feature inclusion counts are pooled over all valid returned alpha/lambda models. The final selected alpha can vary between outer splits. Output tables retain the column `alpha_fixed` for that selected value; it is an output field, while the input configuration uses `alpha_grid`.

With the defaults, each outer training set fits 11 alpha paths with 50 lambda points each: 550 requested coefficient models for feature-frequency counting. Across 100 outer repeats, this requests 55,000 frequency models, produced by 1,100 glmnet path fits. Final training CV separately fits the 11 alpha candidates on the selected TOP_n features using the same fold assignment, then retains one model per outer split. AUC summaries therefore use up to 100 final models. Frequencies use the actual valid returned models, and path failures or truncation are recorded.

### Command-line run and function-only mode

```bash
# Defaults: input.csv, current working directory
Rscript FOREST.R

# Specify input and output directory
Rscript FOREST.R input.csv formal_run

# Show positional argument usage without fitting
Rscript FOREST.R --help
```

The positional arguments override input/output paths. Other settings use the script's `frequency_default_config()` values, which can be edited at the top of the file for command-line use. Running again in the same output directory requires a different directory or an explicit `overwrite = TRUE` setting.

To load functions without immediately fitting:

```r
options(frequency_oof.functions_only = TRUE)
source("FOREST.R", encoding = "UTF-8")
options(frequency_oof.functions_only = FALSE)

frequency_run <- run_frequency_pipeline(list(
  input_file = "input.csv",
  output_dir = "another_run",
  train_fraction = 0.67,
  TOP_n = 10L,
  alpha_grid = seq(0.5, 1, by = 0.05),
  lambda_ratio_min = 0.001,
  n_repeats = 100L,
  n_lambda = 50L,
  pvalue_cutoff = 0.05,
  seed = 123L
))
```

This function call returns the complete result invisibly; use `frequency_run$results_df` and the corresponding table fields.

### Modeling procedure

1. Generate `n_repeats` stratified random outer splits using the specified fraction. Each sample belongs to exactly one of training or validation within a split.
2. With coverage enabled, inspect all planned validation IDs before fitting. If a sample never appears, replace an extra validation occurrence of a same-group sample with the uncovered sample. This preserves each split's per-class counts and the configured number of repeats. The repair uses IDs and planned membership, without consulting expression values or fitted performance; changes are recorded in `split_adjustments.csv`.
3. In each outer training set, test all original features using the selected differential method and retain raw p-values below the threshold. There is no MAD, missing-rate, or preliminary TOP_n filter. If `pvalue_cutoff = NULL`, skip differential testing/filtering.
4. Calculate medians, means, and standard deviations from training data. Remove columns that cannot be imputed or have zero/nonfinite SD. Apply the training transformation to validation data.
5. On the complete preprocessed outer training set, fit one lambda path for each value in `alpha_grid`. Count how often each feature has an absolute coefficient greater than `coef_tol`, pooling across all alpha/lambda models. The frequency denominator uses the actual valid returned models; truncated paths and failures are recorded.
6. Rank by descending inclusion count, break ties by feature name, and select exactly TOP_n positive-frequency features. An insufficient candidate pool is recorded as a failure; TOP_n is not silently reduced.
7. On the selected TOP_n features, compare every alpha using the same training-only CV fold assignment. Within an alpha, select the minimum complete, finite CV deviance and prefer larger lambda on exact ties. Across alpha values, select the smallest CV deviance and prefer larger alpha on exact ties. The chosen CV object's model is fitted on the complete outer training set. Only after this selection, predict the outer validation set and calculate AUC.
8. Repeat and export all predictions, coefficients, AUCs, curves, feature frequencies, splits, and diagnostics.

Final regularization may set some of the TOP_n input coefficients to zero. `N_Used` is the input count; `N_Final_Nonzero` is the final nonzero count.

Every class needs at least five samples overall, at least three training samples, and at least two validation samples. A training fraction of 0.9 therefore needs at least 11 samples per class. For complete planned coverage, each class must also satisfy `n_repeats * validation_count >= class_count`; the script reports the minimum required repeats when this cannot hold.

### Saved files

| File | Contents |
| --- | --- |
| `oof_predictions.csv` | One row per successful validation prediction, including ID, group, repeat, holdout, parameter fields, and `Tumor_Probability` |
| `oof_predictions_summary.csv` | One row per input sample, in input order: planned validation count, successful prediction count, mean/median/SD/min/max probability, and coverage flag |
| `results_df.csv` | All successful training and validation sample predictions; `Set` identifies TRN/VAL |
| `model_df.csv` | Final selected-lambda coefficients for every successful model, including the intercept and zero coefficients; feature coefficients are on the training-standardized scale |
| `auc_df.csv` | One row per outer attempt: training/validation AUC, selected alpha and lambda, feature counts, warnings, and errors |
| `curve_df.csv` | Training and validation FPR/TPR values on `fpr_grid` for each successful model |
| `indexfilterresult.csv` | AUC mean and 2.5%/97.5% empirical quantiles across successful repeats, plus success/failure counts |
| `selected_features.csv`, `feature_frequency.csv` | Per-split selected features and training inclusion-frequency rankings |
| `differential_results.csv` | Training-set differential test results; P is NA when differential testing is disabled |
| `screen_diagnostics.csv`, `final_diagnostics.csv`, `indexfilter_errors.csv` | Frequency-path status; per-alpha training CV deviance, selected candidate flag, warnings/errors; failed outer attempts |
| `split_membership.csv`, `split_adjustments.csv` | Every sample's planned training/validation membership and any coverage-repair swaps |
| `run_config.R`, `session_info.txt`, `run_log.txt`, `run_status.txt` | Configuration, R/package versions, fitting messages, input checksum, and completion/coverage counts |
| `frequency_analysis.rds` | Complete result object and saved models, written when `save_models = TRUE` |

### Interpretation of out-of-fold probabilities

`oof_predictions.csv` contains only **outer validation predictions**. A sample's label and features were excluded from the training, differential analysis, frequency selection, preprocessing estimation, and alpha/lambda selection for the model producing that row. Training-CV `fit.preval` and in-sample training predictions are not used to populate this file.

This is repeated holdout validation, rather than a single K-fold partition. A sample can have multiple valid OOF predictions from different models. `Tumor_Probability_Mean` averages only that sample's successful validation predictions. The SD summarizes their spread and is NA when only one valid prediction exists. Do not pool repeated prediction rows as if they were independent subjects.

For one probability per sample:

```r
scores <- read.csv(
  "oof_predictions_summary.csv",
  colClasses = c(ID = "character"),
  check.names = FALSE
)
scores[, c("ID", "Group", "N_Validation_Predictions", "Tumor_Probability_Mean")]
```

ROC grid interpolation is for curve export. Reported AUCs come from the original predictions using pROC; they are not recalculated from the interpolated grid. The AUC quantile limits summarize repeat variability and are not confidence intervals for a mean. Training CV selects alpha and lambda after training-only feature selection; its deviance is not an independent estimate of the entire selection workflow's performance. Freeze formal-run settings beforehand. When settings were selected using the same cohort's validation results, independent data are needed for a final performance estimate.

### Failures, incomplete coverage, and reruns

Failed outer attempts remain in `auc_df.csv` and `indexfilter_errors.csv`, with missing AUC and an error message. Their predictions are not generated. Constant/unusable features, too few p-significant or positive-frequency features, numerical path failures, or incompatible sample counts can prevent fitting.

Planned validation coverage does not guarantee successful fitting. If a sample never receives a valid validation prediction, it is retained in the summary with `Has_OOF = FALSE`, count zero, and NA probabilities. The default raises an error **after saving the available tables and diagnostics**; `run_status.txt` records `INCOMPLETE_OOF`. No training probability or artificial zero is used to fill it.

Review the errors, then choose a new output directory and correct the relevant configuration. Options include changing the training fraction, increasing repeats when coverage is inadequate, revisiting TOP_n, or explicitly disabling p-value filtering if that is the intended modeling design. `require_complete_oof = FALSE` allows a partial result to return, with the same missing-value and status records.

A complete OOF status means every sample has at least one successful validation prediction. Some outer attempts can still have failed; check the success/failure counts separately. Keep `overwrite = FALSE` for separate runs, or set TRUE when intentionally replacing all named result files.

Saved models can predict a new feature matrix using the same stored training transformation:

```r
saved <- readRDS("frequency_analysis.rds")
model <- saved$analysis$models[[1L]]
new_samples <- read.csv(
  "new_samples.csv", colClasses = c(ID = "character"), check.names = FALSE
)
tumor_probability <- predict_frequency_model(model, new_samples)
```

The CSV reader requires both groups for the training entry point. For unlabeled future samples, provide a numeric data frame with the model's named feature columns directly to `predict_frequency_model()`.

### Modeling documentation and validation status

- [glmnet cross-validation documentation](https://glmnet.stanford.edu/reference/cv.glmnet.html)
- [pROC documentation](https://xrobin.github.io/pROC/)
- [R CSV input documentation](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/read.table.html)

The deliverables were checked for input schema, preservation of the original README, script syntax structure, and planned-validation coverage logic. This revision also checks the complete-training-set alpha/lambda frequency scope, requested-model denominator, training-only final selection, and per-holdout selected-alpha mapping into exported tables. The packaging environment did not contain R/Rscript, so end-to-end glmnet fitting was not executed there. The example input is ready for execution in an R environment with the stated dependencies.


The deliverables were checked for input schema, preservation of the original README, script syntax structure, and planned-validation coverage logic. This revision also checks that the shared repeat setting reaches the outer loop, inner subsample plan, frequency loop, requested-model denominator, configuration, and entry point. The packaging environment did not contain R/Rscript, so end-to-end glmnet fitting was not executed there. The example input is ready for execution in an R environment with the stated dependencies.

