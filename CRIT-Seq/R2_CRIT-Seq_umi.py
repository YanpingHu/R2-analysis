# -*- coding: utf-8 -*-
import re
import gzip
import sys
import subprocess
import os

def F_to_R(sequence):
    transtab = str.maketrans("ACGTacgt", "TGCATGCA")
    return sequence.translate(transtab)[::-1]

def fq(file):
    if re.search('.gz$', file):
        fastq = gzip.open(file, 'rb')
    else:
        fastq = open(file, 'r')
    with fastq as f:
        while True:
            l1 = f.readline().decode()
            if not l1:
                break
            l2 = f.readline().decode()
            l3 = f.readline().decode()
            l4 = f.readline().decode()
            yield [l1, l2, l3, l4]

def genome_umi_qual(fastq, index, R1_start, R2_start):
    index_R = F_to_R(index)
    R1_start_R = F_to_R(R1_start)
    R2_start_R = F_to_R(R2_start)

    fastq_seq = fastq[1].strip()

    index_start = fastq_seq.find(index)
    R1_start_start = fastq_seq.find(R1_start)
    R2_start_start = fastq_seq.find(R2_start)

    index_R_start = fastq_seq.find(index_R)
    R1_start_R_start = fastq_seq.find(R1_start_R)
    R2_start_R_start = fastq_seq.find(R2_start_R)

    if index_start != -1:
        if index_start == R1_start_start -10 - len(index):
            umi_start = index_start + len(index)
            umi_end = R1_start_start
            umi = fastq_seq[umi_start:umi_end]

            start = R1_start_start + len(R1_start)

            if R2_start_start != -1:
                end = R2_start_start
                line = fastq_seq[start:R2_start_start]
            else:
                end = len(fastq_seq)
                line = fastq_seq[start:]
        else:
            umi = ''
            start = ''
            end = ''
            line = ''
        direction = 'R'

    elif R2_start_R_start != -1:
        start = R2_start_R_start + len(R2_start)
        if R1_start_R_start != -1:
            umi_start = R1_start_R_start + len(R1_start)
            umi_end = umi_start + 10
            if umi_end <= len(fastq_seq) and umi_start <= len(fastq_seq):
                umi = F_to_R(fastq_seq[umi_start:umi_end])
            elif umi_end > len(fastq_seq) and umi_start <= len(fastq_seq):
                umi = F_to_R(fastq_seq[umi_start:len(fastq_seq)])
            else:
                umi = ''
            end = R1_start_R_start
            line = fastq_seq[start:end]
        else:
            umi = ''
            end = len(fastq_seq)
            line = fastq_seq[start:]
        direction = 'F'
    else:
        start = ''
        end = ''
        line = ''
        umi = ''
        direction = ''

    if len(line) >= 10:
        line_start = line[0:10]
    else:
        line_start = line


    name = fastq[0].rstrip().split()[0]
    if len(line) > 10:
        qual = fastq[3].strip()[start:end]
    else:
        line = ''
        qual = ''

    return name, umi, line_start, line, qual, direction

def umitag(read1, read2, index, R1_start, R2_start, outfile1, outfile2):
    r1_umitagged = open(outfile1, 'w')
    r2_umitagged = open(outfile2, 'w')
    for r1, r2 in zip(fq(read1), fq(read2)):
        name_1, umi_1, start_1, seq_1, qual_1, direction_1 = genome_umi_qual(r1, index, R1_start, R2_start)
        name_2, umi_2, start_2, seq_2, qual_2, direction_2 = genome_umi_qual(r2, index, R1_start, R2_start)

        if len(umi_1) > 0 or len(umi_2) > 0:
            if len(umi_1) == len(umi_2):
                if umi_1 == umi_2:
                    umi_sum = umi_1
                else:
                    umi_sum = ""
            elif len(umi_1) > len(umi_2):
                umi_sum = umi_1
            else:
                umi_sum = umi_2
        else:
            umi_sum = ''

        if direction_1 == 'F':
            direction_sum_1 = start_1
            direction_sum_2 = start_2
            read1 = seq_1
            read2 = seq_2
            quality1 = qual_1
            quality2 = qual_2
            re_name1 = name_1
            re_name2 = name_2
        else:
            direction_sum_1 = start_2
            direction_sum_2 = start_1
            read1 = seq_2
            read2 = seq_1
            quality1 = qual_2
            quality2 = qual_1
            re_name1 = name_2
            re_name2 = name_1

        if len(umi_sum) == 10 :
            molecular_id = '%s' % (umi_sum)
            read_name_1 = '%s %s' % (re_name1, molecular_id)
            read_name_2 = '%s %s' % (re_name2, molecular_id)

            if len(seq_1) > 0 and len(seq_2) > 0:
                r1_umitagged.write(read_name_1 + '\n' + read1 + '\n' + '+' + '\n' + quality1 + '\n')
                r2_umitagged.write(read_name_2 + '\n' + read2 + '\n' + '+' + '\n' + quality2 + '\n')
    r1_umitagged.close()
    r2_umitagged.close()

def usage():
    print('Usage: python3 script.py read1, read2, index1, R1_start, index2, tmp1, tmp2, read1_out, read2_out')


def main():
    read1 = sys.argv[1]
    read2 = sys.argv[2]
    index = sys.argv[3].upper()
    R1_start = sys.argv[4].upper()
    R2_start = sys.argv[5].upper()
    outfile1 = sys.argv[6]
    outfile2 = sys.argv[7]
    read1_out= sys.argv[8]
    read2_out = sys.argv[9]

    umitag(read1, read2, index, R1_start, R2_start, outfile1, outfile2)

    cmd = 'cat ' + outfile1 + ' | paste - - - - | sort -k2,2 -k1,1 | tr "\t" "\n" >' + read1_out
    subprocess.check_call(cmd, shell=True, env=os.environ.copy())
    cmd = 'cat ' + outfile2 + ' | paste - - - - | sort -k2,2 -k1,1 | tr "\t" "\n" >' + read2_out
    subprocess.check_call(cmd, shell=True, env=os.environ.copy())
    os.remove(outfile1)
    os.remove(outfile2)

try:
    main()
except IndexError:
    usage()

