# 1. Assess raw read quality
fastqc -t 60 input_R1.fastq.gz input_R2.fastq.gz
multiqc .

# 2. Trim adapters and filter low-quality reads
java -jar trimmomatic.jar PE -threads 40 \
    input_R1.fastq.gz input_R2.fastq.gz \
    clean_R1_paired.fastq.gz clean_R1_unpaired.fastq.gz \
    clean_R2_paired.fastq.gz clean_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:5 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:36

# 3. Build the reference genome index (run once per reference)
bowtie2-build --threads 60 reference.fa reference_index

# 4. Align paired-end reads, filter alignments, and sort by coordinates
bowtie2 --threads 20 --end-to-end --very-sensitive --no-unal \
    -x reference_index \
    -1 clean_R1_paired.fastq.gz \
    -2 clean_R2_paired.fastq.gz \
    | samtools view -b -q 10 -F 4 \
    | samtools sort -@ 20 -l 5 -m 10G -o sorted.bam

# 5. Add read group information
picard AddOrReplaceReadGroups \
    I=sorted.bam O=sorted_RG.bam \
    RGID=sample RGLB=library RGPL=ILLUMINA RGPU=unit RGSM=sample

# 6. Remove duplicate reads
picard MarkDuplicates \
    I=sorted_RG.bam O=deduplicated.bam \
    M=duplicate_metrics.txt \
    REMOVE_DUPLICATES=true