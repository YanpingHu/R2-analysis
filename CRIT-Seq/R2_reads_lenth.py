# encoding:utf-8
import re
import gzip
import sys

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

def is_extract(r1,r2,num):
    r1 = r1[1].strip()
    r2 = r2[1].strip()

    if len(r1) >= num and len(r2) >= num:
        return True
    else:
        return False

def plasmid_extract(read1, read2, num, R1_file, R2_file):
    R1_file = open(R1_file, 'w')
    R2_file = open(R2_file, 'w')
    for r1, r2 in zip(fq(read1), fq(read2)):
        result = is_extract(r1,r2,num)
        if result:
            R1_file.write(r1[0].strip() + '\n' + r1[1].strip() + '\n' + '+' + '\n' + r1[3].strip() + '\n')
            R2_file.write(r2[0].strip() + '\n' + r2[1].strip() + '\n' + '+' + '\n' + r2[3].strip() + '\n')

    R1_file.close()
    R2_file.close()

def main():
    read1 = sys.argv[1]
    read2 = sys.argv[2]
    num = int(sys.argv[3])
    R1_file = sys.argv[4]
    R2_file = sys.argv[5]
    plasmid_extract(read1, read2, num, R1_file, R2_file)

if __name__ == '__main__':
    main()
