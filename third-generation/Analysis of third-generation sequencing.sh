# 1. The resulting reads.bam files containing HiFi Reads were converted to reads.fastq format using samtools software with bam2fq parameters.
id="name_of_bam"
samtools bam2fq ${id}.bam > ${id}.fastq 

# 2. Then, reads that could be mapped onto the sequences of 200 bp LHA-RHA or 3' UTR were extracted, using BWA with the parameters: ‘-x pacbio’ and ‘-T 70’.
index_of_first_mapping='third-generation-rDNA-3-UTR-mapping'
bwa mem -t 20 -x pacbio -T 70 ${index_of_first_mapping} ${id}.fastq -o ${id}.sam 
samtools sort -@ 10 -O bam ${id}.sam -o ${id}.sorted.bam
samtools view -bF 4 ${id}.sorted.bam >${id}-align.bam
samtools view ${id}-align.bam >${id}-align.sam
awk '{print $1}' ${id}-align.sam | sort | uniq >${id}-align-name.txt
seqkit grep -f ${id}-align-name.txt ${id}.fastq > ${id}-align.fastq

# 3. Extracted HiFi reads that aligned with plasmids sequences were discarded.
index_of_plasmid='third-generation-plasmid'
bwa mem -t 20 -x pacbio ${index_of_plasmid} ${id}-align.fastq -o ${id}-plasmid.sam
samtools sort -@ 10 -O bam ${id}-plasmid.sam -o ${id}-plasmid.sorted.bam
samtools view -bF 4 ${id}-plasmid.sorted.bam >${id}-plasmid-align.bam
samtools view ${id}-plasmid-align.bam >${id}-plasmid-align.sam
awk '{print $1}' ${id}-plasmid-align.sam >${id}-plasmid-align-name.txt
python two_file_inter.py ${id}-align-name.txt ${id}-plasmid-align-name.txt ${id}-plasmid-noalign-name.txt
seqkit grep -f ${id}-plasmid-noalign-name.txt ${id}-align.fastq > ${id}-plasmid-noalign.fastq

# 4. The remaining reads, which aligned with the sequence of the 200 bp LHA+RHA with > 100 matches, were considered unperturbed rDNA loci. Reads containing the 3' UTR sequences were regarded as transgene-integrated reads (TIRs).  
python R2_HIFI_positive_align.py sam_extract ${id}-align.sam ${id}-plasmid-noalign-name.txt ${id}-plasmid-noalign-name.sam
python R2_HIFI_positive_align.py hifi_sam_classification ${id}-plasmid-noalign-name.sam LHA-RHA 100 ${id}-negetive-name.txt ${id}-unknown-name.txt ${id}-positive-name.txt
seqkit grep -f ${id}-positive-name.txt ${id}-align.fastq > ${id}-positive.fastq

# 5. Reads containing the ddPCR-amplified GAPDH sequence were used to represent the number of sequenced genomes.
index_GAPDH='third-generation-GAPDH'
bwa mem -t 20 -x pacbio -T 70 ${index_GAPDH} ${id}.fastq  -o ${id}-GAPDH.sam 
samtools sort -@ 10 -O bam ${id}-GAPDH.sam -o ${id}-GAPDH.sorted.bam
samtools view -bF 4 ${id}-GAPDH.sorted.bam >${id}-GAPDH-align.bam
samtools view ${id}-GAPDH-align.bam >${id}-GAPDH-align.sam
awk '{print $1}' ${id}-GAPDH-align.sam  | sort | uniq >${id}-GAPDH-align-name.txt