# encoding:utf-8

import sys
import re
from collections import Counter
def sam_to_dic(sam_file):
    dic = {}

    with open(sam_file, 'r') as sams:
        for sam in sams:
            if not sam.startswith('@'):
                full_read_name, sam_flag, chromosome, position, mapq, cigar, name_of_mate, position_of_mate, template_length, read_sequence, read_quality = sam.strip().split()[
                                                                                                                                                            :11]
                position = int(position)
                if chromosome not in dic.keys():
                    dic[chromosome] = {}
                    dic[chromosome][position] = [cigar]
                else:
                    if position not in dic[chromosome].keys():
                        dic[chromosome][position] = [cigar]
                    else:
                        dic[chromosome][position].append(cigar)
    return dic

def sam_to_dic_name(sam_file):
    dic = {}

    with open(sam_file, 'r') as sams:
        for sam in sams:
            if not sam.startswith('@'):
                full_read_name, sam_flag, chromosome, position, mapq, cigar, name_of_mate, position_of_mate, template_length, read_sequence, read_quality = sam.strip().split()[
                                                                                                                                                            :11]
                position = int(position)
                inf = full_read_name + '-' + chromosome + '-' + str(position)
                if chromosome not in dic.keys():
                    dic[chromosome] = {}
                    dic[chromosome][position] = [inf]
                else:
                    if position not in dic[chromosome].keys():
                        dic[chromosome][position] = [inf]
                    else:
                        dic[chromosome][position].append(inf)
    return dic

def identified_to_dic(identified_file):
    dic = {}

    with open(identified_file, 'r') as identifieds:
        for identified in identifieds:
            if not identified.startswith('#'):
                Chromosome, Position = identified.strip().split()[6:8]
                Position = int(Position)
                if Chromosome not in dic.keys():
                    dic[Chromosome] = {}
                    dic[Chromosome][Position] = identified.strip()
                else:
                    if Position not in dic[Chromosome].keys():
                        dic[Chromosome][Position] = identified.strip()
                    else:
                        dic[Chromosome][Position].append(identified.strip())
    return dic

def identified_add_align_inf(identified_file, sam_file, out_file):
    out_file = open(out_file, 'w')
    dic_id = identified_to_dic(identified_file)
    dic_sam = sam_to_dic(sam_file)
    sum_dic = {}

    for i in dic_id.keys():
        sum_dic[i] = {}
        for j in dic_id[i].keys():
            sum_dic[i][j] = []
            most_frequent_position = int(j)
            min_position = most_frequent_position - 35
            max_position = most_frequent_position + 35
            for ij in dic_sam[i].keys():
                if ij >= min_position and ij <= max_position:
                    sum_dic[i][j].extend(dic_sam[i][ij])
    out_file.write('\t'.join(['#BED Chromosome', 'BED Min.Position',
                       'BED Max.Position', 'BED Name', 'Filename', 'WindowIndex', 'Chromosome', 'Position',
                       'Sequence', '+.mi', '-.mi', 'bi.sum.mi', 'bi.geometric_mean.mi', '+.total',
                       '-.total', 'total.sum', 'total.geometric_mean', 'primer1.mi', 'primer2.mi',
                       'primer.geometric_mean',
                       'position.stdev', 'Off-Target Sequence', 'Mismatches', 'Length', 'BED off-target Chromosome',
                       'BED off-target start', 'BED off-target end', 'BED off-target name', 'BED Score', 'Strand',
                       'Cells', 'Targetsite', 'Target Sequence', 'Cigar']) + '\n')
    for i in dic_id.keys():
        for j in dic_id[i].keys():
            element_counts = Counter(sum_dic[i][j])

            element_count_dict = dict(element_counts)
            sum_cigar = ''

            for cigar in element_count_dict.keys():
                sum_cigar += cigar + "_" + str(element_count_dict[cigar]) + ':'
            out_file.write(dic_id[i][j] + '\t' + sum_cigar + '\n')

    out_file.close()


def identified_add_ID_inf(identified_file, sam_file, out_file):
    out_file = open(out_file, 'w')
    dic_id = identified_to_dic(identified_file)
    dic_sam = sam_to_dic_name(sam_file)
    sum_dic = {}

    for i in dic_id.keys():
        sum_dic[i] = {}
        for j in dic_id[i].keys():
            sum_dic[i][j] = []
            most_frequent_position = int(j)
            min_position = most_frequent_position - 25
            max_position = most_frequent_position + 25
            for ij in dic_sam[i].keys():
                if ij >= min_position and ij <= max_position:
                    sum_dic[i][j].extend(dic_sam[i][ij])
    out_file.write('\t'.join(['#BED Chromosome', 'BED Min.Position',
                       'BED Max.Position', 'BED Name', 'Filename', 'WindowIndex', 'Chromosome', 'Position',
                       'Sequence', '+.mi', '-.mi', 'bi.sum.mi', 'bi.geometric_mean.mi', '+.total',
                       '-.total', 'total.sum', 'total.geometric_mean', 'primer1.mi', 'primer2.mi',
                       'primer.geometric_mean',
                       'position.stdev', 'Off-Target Sequence', 'Mismatches', 'Length', 'BED off-target Chromosome',
                       'BED off-target start', 'BED off-target end', 'BED off-target name', 'BED Score', 'Strand',
                       'Cells', 'Targetsite', 'Target Sequence', 'Full_name']) + '\n')
    for i in dic_id.keys():
        for j in dic_id[i].keys():

            # element_counts = Counter(sum_dic[i][j])
            #

            # element_count_dict = dict(element_counts)
            sum_name = ''

            for name in sum_dic[i][j]:
                sum_name += name  + ':'
            sum_name.strip(':')
            out_file.write(dic_id[i][j] + '\t' + sum_name + '\n')

    out_file.close()

def identified_sumcount_ratio(file, identified_outfile, sumcount, ratio, seq_file, bed_file):
    limit_conut = sumcount * ratio
    identified_outfile = open(identified_outfile, 'w')
    seq_file = open(seq_file, 'w')
    seq_dic = {}
    chr_dic = {}
    idn_dic = {}
    name_dic = {}
    bed_file = open(bed_file, 'w')

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
                all_name = line.strip().split('(')[-1].split(')')[0].lstrip("'").rstrip("'").rstrip("',")
                name_lis = []
                if ':' not in all_name and ',' not in all_name:
                    name_lis = [all_name.strip().strip("'")]
                elif ':' in all_name and ',' not in all_name:
                    hh = all_name.split(':')
                    for mm in hh:
                        name_lis.append(mm.strip().strip("'"))
                elif ':' not in all_name and ',' in all_name:
                    hh = all_name.split(',')
                    for mm in hh:
                        name_lis.append(mm.strip().strip("'"))
                else:
                    hh = all_name.split(':')
                    for mm in hh:
                        ss = mm.split(',')
                        for xx in ss:
                            name_lis.append(xx.strip().strip("'"))

                # print(all_name,name_lis)
                len_lis = len(name_lis)
                if sum_count == len_lis:

                    print(Chromosome, Position, sum_count, len_lis, 'True')
                else:

                    print(Chromosome, Position, sum_count, len_lis, 'False', all_name, name_lis)
                # print(all_name, name_lis)
                # name_lis = line.strip().rstrip(':').split()[-1].split(':')
                #
                #
                # if sum_count >= limit_conut:
                #     identified_outfile.write(line)
                if seq not in seq_dic.keys():
                    chr_dic[seq] = Chromosome + '\t' + str(start) + '\t' + str(end)
                    seq_dic[seq] = sum_count
                    idn_dic[seq] = line
                    name_dic[seq] = []
                    name_dic[seq].extend(name_lis)
                else:
                    seq_dic[seq] += sum_count
                    name_dic[seq].extend(name_lis)
                    # bed_file.write(Chromosome + '\t' + str(start) + '\t' + str(end) + '\t' + str(sum_count) + '\n')
    seq_dic_sort = dict(sorted(seq_dic.items(), key=lambda x: x[1], reverse=True))
    # print(seq_dic_sort)
    for i in seq_dic_sort.keys():
        # print(seq_dic_sort[i])
        if seq_dic_sort[i] >= limit_conut:
            # print(i, name_dic[i])
            seq_file.write('>' + str(seq_dic_sort[i]) + '\n' + i + '\n')
            bed_file.write(chr_dic[i] + '\t' + str(seq_dic_sort[i]) + '\n')
            for j in name_dic[i]:
                identified_outfile.write(chr_dic[i] + '\t' + str(seq_dic_sort[i]) + '\t' + j.strip("'") + '\n')


    identified_outfile.close()
    seq_file.close()
    bed_file.close()

def target_align_inf_extract(id_file, sam_file, out_file):
    out_file = open(out_file, 'w')
    id_dic = {}

    with open(id_file, 'r') as ids:
        for id in ids:
            name = id.strip().split()[-1]
            id_dic[name] = id.strip()

    with open(sam_file, 'r') as sams:
        for sam in sams:
            if not sam.startswith('@'):
                full_read_name, sam_flag, chromosome, position, mapq, cigar, name_of_mate, position_of_mate, template_length, read_sequence, read_quality = sam.strip().split()[
                                                                                                                                                            :11]
                if full_read_name in id_dic.keys():
                    out_file.write(id_dic[full_read_name] + '\t' + cigar + '\n')
    out_file.close()

def is_rDNA(rDNA_inf_file, all_file, rDNA_file, genome_file):
    rDNA_file = open(rDNA_file, 'w')
    genome_file = open(genome_file, 'w')
    rDNA_lis = []

    with open(rDNA_inf_file, 'r') as lls:
        for ll in lls:
            chr_name = ll.strip().split()[0]
            count = ll.strip().split()[1]
            rDNA_lis.append([chr_name, count])

    with open(all_file, 'r') as mms:
        for mm in mms:
            chr_name = mm.strip().split()[0]
            count = mm.strip().split()[3]
            rDNA_result = False
            for i in rDNA_lis:
                if chr_name == i[0] and count == i[1]:
                    rDNA_result = True

            if rDNA_result:
                rDNA_file.write(mm)
            else:
                genome_file.write(mm)


    rDNA_file.close()
    genome_file.close()

def is_MS(s):

    numbers = re.findall(r'\d+', s)
    letters = re.findall(r'[a-zA-Z]+', s)
    # print(numbers)
    # print(letters)
    if letters == ['M', 'S']:
        return int(numbers[0]), int(numbers[0])  + int(numbers[1])
    else:
        return False, False

def sam_split_S_1(sam_file, outfile, namefile):
    name_lis = {}
    outfile = open(outfile, 'w')

    with open(namefile, 'r') as lls:
        for ll in lls:
            name_lis[ll.strip().split()[4]] = ll.strip()

    with open(sam_file, 'r') as sams:
        for sam in sams:
            if not sam.startswith('@'):
                full_read_name, sam_flag, chromosome, position, mapq, cigar, name_of_mate, position_of_mate, template_length, read_sequence, read_quality = sam.strip().split()[
                                                                                                                                                            :11]

                if full_read_name in name_lis.keys():
                    if chromosome == '*':
                        pass
                    else:
                        start, end = is_MS(cigar)
                        if start or end:
                            seq_S = read_sequence[start:end][0]
                            outfile.write(name_lis[full_read_name] + '\t' + seq_S + '\n')
    outfile.close()

def split_string_to_dict(s):

    numbers = re.findall(r'\d+', s)
    letters = re.findall(r'[a-zA-Z]+', s)
    dic = {}

    for i in range(len(letters)):
        if letters[i] not in dic.keys():
            dic[letters[i]] = []
            dic[letters[i]].append(int(numbers[i]))
        else:
            dic[letters[i]].append(int(numbers[i]))

    return dic
def sam_extract(sam_file, name_file, outfile):
    name_lis = []
    outfile = open(outfile, 'w')
    dic = {}

    with open(name_file, 'r') as lls:
        for ll in lls:
            name_lis.append(ll.strip().split()[4])

    with open(sam_file, 'r') as sams:
        for sam in sams:
            if sam.startswith('@'):
                outfile.write(sam)
            else:
                full_read_name, sam_flag, chromosome, position, mapq, cigar, name_of_mate, position_of_mate, template_length, read_sequence, read_quality = sam.strip().split()[:11]
                cigar_dic = split_string_to_dict(cigar)
                if 'M' in cigar_dic.keys():
                    sum_match = sum(cigar_dic['M'])
                else:
                    sum_match = 0
                if full_read_name in name_lis:
                    if full_read_name not in dic.keys():
                        dic[full_read_name] = [sum_match, sam.strip()]
                    else:
                        if sum_match > dic[full_read_name][0]:
                            dic[full_read_name] = [sum_match, sam.strip()]
    print(len(list(dic.keys())))
    for i in dic.keys():
        outfile.write(dic[i][1] + '\n')
    outfile.close()

def main():

    parament = sys.argv[1]
    if parament == 'identified_add_align_inf':
        identified_file, sam_file, out_file = sys.argv[2:5]
        identified_add_align_inf(identified_file, sam_file, out_file)
    elif parament == 'identified_add_ID_inf':
        identified_file, sam_file, out_file = sys.argv[2:5]
        identified_add_ID_inf(identified_file, sam_file, out_file)
    elif parament == 'identified_sumcount_ratio':
        file, identified_outfile, sumcount, ratio, seq_file, bed_file = sys.argv[2:8]
        sumcount = int(sumcount)
        ratio = float(ratio)
        identified_sumcount_ratio(file, identified_outfile, sumcount, ratio, seq_file, bed_file)
    elif parament == 'target_align_inf_extract':
        id_file, sam_file, out_file = sys.argv[2:5]
        target_align_inf_extract(id_file, sam_file, out_file)
    elif parament == 'is_rDNA':
        rDNA_inf_file, all_file, rDNA_file, genome_file = sys.argv[2:6]
        is_rDNA(rDNA_inf_file, all_file, rDNA_file, genome_file)
    elif parament == 'sam_split_S_1':
        sam_file, outfile, namefile = sys.argv[2:5]
        sam_split_S_1(sam_file, outfile, namefile)
    elif parament == 'sam_extract':
        sam_file, name_file, outfile = sys.argv[2:5]
        sam_extract(sam_file, name_file, outfile)
if __name__ == '__main__':
    main()
