# 1. Quality control and filtering
fastp -w 16 \
    --qualified_quality_phred 30 \
    --length_required 36 \
    --cut_mean_quality 30 --cut_front --cut_tail \
    --trim_front1 5 --trim_tail1 3 \
    --trim_front2 5 --trim_tail2 3 \
    -i input_R1.fq.gz -I input_R2.fq.gz \
    -o clean_R1.fq.gz -O clean_R2.fq.gz \
    -h fastp.html -j fastp.json

# 2. Build the reference genome index (run once per reference)
bwa-mem2 index reference.fa

# 3. Align paired-end reads
bwa-mem2 mem -t 16 reference.fa clean_R1.fq.gz clean_R2.fq.gz > aligned.sam

# 4. Sort alignments by genomic coordinates
samtools sort -@ 20 -l 6 -m 10G aligned.sam -o sorted.bam

# 5. Add read group information
picard AddOrReplaceReadGroups \
    I=sorted.bam O=sorted_RG.bam \
    RGID=sample RGLB=library RGPL=ILLUMINA RGPU=unit RGSM=sample

# 6. Remove duplicate reads
gatk MarkDuplicates \
    I=sorted_RG.bam O=deduplicated.bam \
    M=duplicate_metrics.txt \
    REMOVE_DUPLICATES=true
    