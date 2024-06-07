# 1. To evaluate the on-target indel ratio, primer sequences were first removed from reads. 
id='name_of_amplicon-seq'
python R2_amplicon-seq.py ${id}_R1.fq.gz ${id}_R2.fq.gz CAAAGCATCGCGAAGGCCCG TGGCGGAATCAGCGGGGAAA ${id}_R1.fq ${id}_R2.fq ${id}_R1-tmp ${id}_R2-tmp &

# 2. Then, CRISPResso2 software was used to quantitatively analyze indels, with the parameters: '--quantification_window_size 1' and '--ignore_substitutions'. 
CRISPResso --fastq_r1 ${id}_R1.fq --fastq_r2 ${id}_R2.fq --ignore_substitutions --guide_seq ctatgactctcttaaggtag --amplicon_seq cggcgggtgttgacgcgatgtgatttctgcccagtgctctgaatgtcaaagtgaagaaattcaatgaagcgcgggtaaacggcgggagtaactatgactctcttaaggtagccaaatgcctcgtcatctaattagtgacgcgcatgaatggatgaacgagattcccactgtccctacctactatccagcgaaaccacagccaagggaacgggct --quantification_window_size 1 -n ${id}_withumi_1_substitutions -o CRISPResso &
