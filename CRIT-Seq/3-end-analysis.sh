# 1. Sequenced reads were split based on barcode and 3’ UTR sequences.
name='name_of_sequencing_read'
file='split_file_name'
python R2_CRIT-Seq_umi.py ${name}_R1.fq.gz ${name}_R2.fq.gz CACTAG ATGCCGTTCTCGCTACT ATCGAACACACATAAATTATATATGTTCCTTTTCCGGGTTAAGTAAAGGTGGCC ${file}-r1 ${file}-r2 ${file}-R1.fq ${file}-R2.fq 
gzip ${file}-R1.fq ${file}-R2.fq 
python R2_reads_lenth.py ${file}-R1.fq.gz ${file}-R2.fq.gz 103 ${file}-len_R1.fq ${file}-len_R2.fq
gzip ${file}-len_R1.fq ${file}-len_R2.fq

# 2. Reads aligned with plasmid sequences used in this study were discarded to exclude plasmids contamination.
index_plasmid='3-end-plasmid'
bwa mem -t 20 -T 15 -k 15 ${index_plasmid} ${file}-len_R1.fq.gz ${file}-len_R2.fq.gz -o ${file}-plasmid.sam 
samtools sort -@ 10 -O bam ${file}-plasmid.sam -o ${file}-plasmid.sorted.bam
samtools view -bF 4 ${file}-plasmid.sorted.bam >${file}-plasmid-align.bam
samtools view ${file}-plasmid-align.bam >${file}-plasmid-align.sam
awk '{print $1}' ${file}-plasmid-align.sam | sort | uniq >${file}-plasmid-align-name.txt
cat ${file}-plasmid-align-name.txt | while read id
do
	echo "${id}/1" >>${file}-plasmid-align-name_1.txt 
	echo "${id}/2" >>${file}-plasmid-align-name_2.txt 
done
zcat ${file}-len_R1.fq.gz | awk 'NR%4==1'  | awk '{print $1}' | awk -F'@' '{print $2}' >${file}-align-name_1.txt
zcat ${file}-len_R2.fq.gz | awk 'NR%4==1'  | awk '{print $1}' | awk -F'@' '{print $2}' >${file}-align-name_2.txt
python two_file_inter.py ${file}-align-name_1.txt ${file}-plasmid-align-name_1.txt ${file}-plasmid-noalign-name_1.txt
python two_file_inter.py ${file}-align-name_2.txt ${file}-plasmid-align-name_2.txt ${file}-plasmid-noalign-name_2.txt
seqkit grep -f ${file}-plasmid-noalign-name_1.txt ${file}-len_R1.fq.gz  > ${file}-plasmid-noalign_R1.fastq
seqkit grep -f ${file}-plasmid-noalign-name_2.txt ${file}-len_R2.fq.gz  > ${file}-plasmid-noalign_R2.fastq

# 3. Then, PCR-generated duplicate reads were consolidated based on unique molecular identifier sequences. 
python R2_consolidate.py ${file}-plasmid-noalign_R1.fastq ${file}_R1_out.consolidated.fastq 14 0.8 
python R2_consolidate.py ${file}-plasmid-noalign_R2.fastq ${file}_R2_out.consolidated.fastq 14 0.8 

# 4. For on-target read mapping, sequences comprising residual 3’ UTR, 100 bp RHA, and extra 1,000 bp DNA downstream of RHA in 219 28s rDNA genes, served as on-target reference sequences (ORSs). Consolidated reads shorter than the length of the residual 3’ UTR plus RHA were discarded, and the remaining reads aligned to the ORSs using BWA software. Only reads that aligned across the 100 bp RHA sequence were defined as on-target. 
index_positive='3-end-on-target-reference-sequences'
bwa mem -t 20 -M ${index_positive} ${file}_R1_out.consolidated.fastq 1>${file}-positive.sam 2>${file}-positive.align.log 
samtools sort -@ 10 -O bam ${file}-positive.sam -o ${file}-positive.sorted.bam
samtools view -bF 4 ${file}-positive.sorted.bam >${file}-positive-align.bam
samtools view ${file}-positive-align.bam >${file}-positive-align.sam
python R2_CRIT-Seq_positive_align.py guideseq_sam_positive ${file}-positive.sam ${file}-positive-result.txt 103
awk '{print $1}' ${file}-positive-result.txt >${file}-positive-name.txt

# 5. The ORS-unmatched portion in each remaining read was mapped onto the Telomere-to-Telomere (T2T) human genome using BWA to identify peaks. any reads that could not be mapped to the T2T genome were discarded.
index_t2t='index_of_t2t'
genome_t2t='genome_of_t2t'
awk 'NR%4==1' ${file}_R1_out.consolidated.fastq | awk '{print $1}' | awk -F'@' '{print $2}' >${file}-all-name.txt
python two_file_inter.py ${file}-all-name.txt ${file}-positive-name.txt ${file}-negative-name.txt
python R2_CRIT-Seq_positive_align.py sam_split_S ${file}-positive.sam ${file}-negative_R1.fq ${file}-negative-name.txt
bwa mem -t 20 -M ${index_t2t} ${file}-negative_R1.fq 1>${file}-genome.sam 2>${file}-genome.align.log 
samtools sort -@ 10 -O bam ${file}-genome.sam -o ${file}-genome.sorted.bam
samtools view -bF 4 ${file}-genome.sorted.bam >${file}-genome-align.bam
samtools view ${file}-genome-align.bam >${file}-genome-align.sam
python R2_CRIT-Seq_positive_align.py choose_max_matches ${file}-genome-align.sam ${file}-genome-max-matches-align.sam
python R2_CRIT-Seq_peak.py ${file}-genome-max-matches-align.sam ${genome_t2t} ./identified/${file}_peak.txt ${file} ${file} ''

# 6. A peak, only when its contained reads number was over 0.1% of total sequenced reads (> 1,000×), was retained and defined as a potential off-target region (POR). 
count='reads_number_of_3end'
# ${file}_peak_0.001.bed uses to annotate using the ChIPseeker package
python R2_CRIT-Seq_off_target.py identified_sumcount_ratio ./identified/${file}_peak.txt ${file}_peak_0.001.txt ${count} 0.001 ${file}_peak_0.001.fasta ${file}_peak_0.001.bed
python R2_CRIT-Seq_off_target.py to_simple_file ${file}_peak_0.001.fasta ${file} GGTAAACGGCGGGAGTAACTATGACTCTCTTAAGGTAGCCAAATGCCTCGTCATCTAATTAGTGACGCGC
ls ${file}*.fasta | awk -F'.fasta' '{print $1}' | while read id
do
	mafft --maxiterate 1000 --localpair  ${id}.fasta > ${id}.fasta 
done

# 7. To analyze the precision of 3'-end integration at28s rDNA loci, the indels occurring at 3’ junction were calculated using CRISPResso2 software with the parameters: “quantification_window_size 1 -amas 0 --default_min_aln_score 0 --exclude_bp_from_right 0 --plot_window_size 3”, by inputting on-target reads mapping to ORSs. 
seqkit grep -f ${file}-positive-name.txt ${file}_R1_out.consolidated.fastq > ${file}-positive_R1.fastq
CRISPResso --fastq_r1 ${file}-positive_R1.fastq --ignore_substitutions --guide_seq gacgaggcatttggctattt --amplicon_seq GTGAAAAAACGATCAGAGTAGTGGTATTTCACCGGCGGCCCGCAGGGCCGGCGGACCCCGCCCCGGGCCCCTCGCGGGGACACCGGGGGGGCGCCGGGGGCCTCCCACTTATTCTACACCTCTCATGTCTCTTCACCGTGCCAGACTAGAGTCAAGCTCAACAGGGTCTTCTTTCCCCGCTGATTCCGCCAAGCCCGTTCCCTTGGCTGTGGTTTCGCTGGATAGTAGGTAGGGACAGTGGGAATCTCGTTCATCCATTCATGCGCGTCACTAATTAGATGACGAGGCATTTGGCTATTT --quantification_window_size 1 -amas 0 --default_min_aln_score 0 --exclude_bp_from_right 0 --plot_window_size 3 -n ${file}_substitutions -o CRISPResso

