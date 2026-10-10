# 1. Assess raw read quality
fastqc -t 60 input_R1.fq.gz input_R2.fq.gz
multiqc .

# 2. Trim adapters and filter low-quality reads
java -jar trimmomatic.jar PE -threads 40 \
    input_R1.fq.gz input_R2.fq.gz \
    clean_R1_paired.fastq.gz clean_R1_unpaired.fastq.gz \
    clean_R2_paired.fastq.gz clean_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:5 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:36

# 3. Build the STAR genome index (run once per reference)
STAR --runMode genomeGenerate \
    --runThreadN 40 \
    --genomeDir genome_index \
    --genomeFastaFiles reference.fa \
    --sjdbGTFfile annotation.gtf

# 4. Align paired-end RNA-seq reads
STAR --runMode alignReads \
    --runThreadN 40 \
    --readFilesCommand zcat \
    --genomeDir genome_index \
    --readFilesIn clean_R1_paired.fastq.gz clean_R2_paired.fastq.gz \
    --quantMode TranscriptomeSAM GeneCounts \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix sample_

# 5. Sort genomic alignments by coordinates
# Use genomic alignments for featureCounts with a genomic GTF annotation
samtools sort -@ 20 -l 6 -m 10G \
    sample_Aligned.out.bam \
    -o sorted.bam

# 6. Sort transcriptome alignments if needed
samtools sort -@ 20 -l 6 -m 10G \
    sample_Aligned.toTranscriptome.out.bam \
    -o transcriptome_sorted.bam

# 7. Add read group information to genomic alignments
picard AddOrReplaceReadGroups \
    I=sorted.bam O=sorted_RG.bam \
    RGID=sample RGLB=library RGPL=ILLUMINA RGPU=unit RGSM=sample

# 8. Remove duplicate reads
picard MarkDuplicates \
    I=sorted_RG.bam O=deduplicated.bam \
    M=duplicate_metrics.txt \
    REMOVE_DUPLICATES=true

# 9. Count exon-overlapping alignments by gene
featureCounts -T 32 -t exon -g gene_id -p \
    -a annotation.gtf -F GTF \
    -o gene_counts.txt \
    deduplicated.bam

# 10. Count exon-overlapping alignments by transcript annotation
featureCounts -T 36 -t exon -g transcript_id -p \
    -a annotation.gtf -F GTF \
    -o transcript_counts.txt \
    deduplicated.bam

# 11. Remove the metadata comment from the transcript count table
sed -i '1d' transcript_counts.txt