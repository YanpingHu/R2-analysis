import sys
import pandas as pd

def merge_htseq_count_genename(file,outfile):
    dic = {}
    outfile = open(outfile, 'w')
    outfile.write(',')

    with open(file, 'r') as lines:
        for line in lines:
            fpkm_file = line.strip()
            all_a = False
            with open(fpkm_file, 'r') as lls:
                for ll in lls:
                    if ll.startswith('A1BG'):
                        all_a = True

                    if all_a:
                        gene_name = ll.strip().split()[0]
                        num = ll.strip().split()[1]
                        if gene_name not in dic.keys():
                            dic[gene_name] = []
                            dic[gene_name].append(num)
                        else:
                            dic[gene_name].append(num)
            outfile.write(fpkm_file + ',')
    outfile.write('\n')
    for i in dic.keys():
        if not i.startswith('__'):
            outfile.write(i + ',' )
            for j in dic[i]:
                outfile.write(str(j) + ',')
            outfile.write('\n')
    outfile.close()
# merge_htseq_count('name.txt','htseq_count_result.txt')

def main():
    file, outfile = sys.argv[1:3]
    merge_htseq_count_genename(file, outfile)

if __name__ == '__main__':
    main()


