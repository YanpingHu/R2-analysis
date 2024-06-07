# encoding:utf-8
import sys

def identified_to_bed(file, outfile):
    outfile = open(outfile, 'w')

    with open(file, 'r') as lines:
        for line in lines:
            if line.startswith('#'):
                pass
            else:
                Chromosome = line.strip().split()[6]
                Position = int(line.strip().split()[7])
                sum_count = line.strip().split()[11]
                start = Position - 35
                end = Position + 35
                outfile.write(Chromosome + '\t' + str(start) + '\t' + str(end) + '\t' + str(sum_count) + '\n')

    outfile.close()

def identified_sumcount_ratio(file, identified_outfile, sumcount, ratio, seq_file, bed_file):
    limit_conut = sumcount * ratio
    identified_outfile = open(identified_outfile, 'w')
    seq_file = open(seq_file, 'w')
    seq_dic = {}
    chr_dic = {}
    idn_dic = {}
    bed_file = open(bed_file, 'w')

    identified_outfile.write('\t'.join(['#BED Chromosome', 'BED Min.Position',
                           'BED Max.Position', 'BED Name', 'Filename', 'WindowIndex', 'Chromosome', 'Position',
                           'Sequence', '+.mi', '-.mi', 'bi.sum.mi', 'bi.geometric_mean.mi', '+.total',
                           '-.total', 'total.sum', 'total.geometric_mean', 'primer1.mi', 'primer2.mi',
                           'primer.geometric_mean',
                           'position.stdev', 'Off-Target Sequence', 'Mismatches', 'Length', 'BED off-target Chromosome',
                           'BED off-target start', 'BED off-target end', 'BED off-target name', 'BED Score', 'Strand',
                           'Cells', 'Targetsite', 'Target Sequence']) + '\n')

    with open(file, 'r') as lines:
        for line in lines:
            if line.startswith('#'):
                pass
            else:
                Chromosome = line.strip().split()[6]
                Position = int(line.strip().split()[7])
                sum_count = int(line.strip().split()[11])
                seq = line.strip().split()[8]
                start = Position - 35
                end = Position + 35
                #
                #
                # if sum_count >= limit_conut:
                #     identified_outfile.write(line)
                if seq not in seq_dic.keys():
                    chr_dic[seq] = Chromosome + '\t' + str(start) + '\t' + str(end)
                    seq_dic[seq] = sum_count
                    idn_dic[seq] = line
                else:
                    seq_dic[seq] += sum_count
                    # bed_file.write(Chromosome + '\t' + str(start) + '\t' + str(end) + '\t' + str(sum_count) + '\n')
    seq_dic_sort = dict(sorted(seq_dic.items(), key=lambda x: x[1], reverse=True))
    # print(seq_dic_sort)
    for i in seq_dic_sort.keys():
        # print(seq_dic_sort[i])
        if seq_dic_sort[i] >= limit_conut:
            seq_file.write('>' + str(seq_dic_sort[i]) + '\n' + i + '\n')
            bed_file.write(chr_dic[i] + '\t' + str(seq_dic_sort[i]) + '\n')
            identified_outfile.write(idn_dic[i])

    identified_outfile.close()
    seq_file.close()
    bed_file.close()




def main():
    parameter = sys.argv[1]
    if parameter == 'identified_to_bed':
        file, outfile = sys.argv[2:4]
        identified_to_bed(file, outfile)

    elif parameter == 'identified_sumcount_ratio':
        file, identified_outfile, sumcount, ratio, seq_file, bed_file = sys.argv[2:8]
        sumcount = int(sumcount)
        ratio = float(ratio)
        identified_sumcount_ratio(file, identified_outfile, sumcount, ratio, seq_file, bed_file)

if __name__ == '__main__':
    main()
