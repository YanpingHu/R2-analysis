# 1. Reads from RNA-seq experiments were aligned with the T2T human genome using Hisat2 (v2.2.1) software.
t2t='t2t-hisat-index'
sample='name_of_RNA-seq_sample'
hisat2 -p 18 --dta-cufflinks -q -x ${t2t} -1 ${sample}_R1.fq.gz -2 ${sample}_R2.fq.gz -S ${sample}.sam
samtools view -Su -q 30 ${sample}.sam | samtools sort -@ 18 - > ${sample}.sorted.bam
samtools index ${sample}.sorted.bam

# 2. Read counts were calculated using htseq-count. 
gtf='t2t_gtf_file'
htseq-count -f bam ${sample}.sorted.bam ${gtf} > ${sample}.HTSeq.out
ls *.HTSeq.out >inf.txt

python htseq_count_merge.py inf.txt sum_count_result.csv
sed -i 's@,$@@g' sum_count_result.csv 


# 3. By comparing the read counts of experimental and control groups using the DEseq2 package, significant differentially expressed genes were defined as having a FDR-adjusted p value < 0.05 and a log2 fold change > 1. Volcano plots were generated in R 4.3.1 using ggplot2 package.
DEseq2_Volcano.R
