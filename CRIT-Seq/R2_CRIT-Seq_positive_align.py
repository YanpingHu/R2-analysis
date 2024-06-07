# encoding:utf-8

import re
import sys

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

def guideseq_sam_positive(sam_file, outfile, matches_number):
    outfile = open(outfile, 'w')
    dic = {}
    matches_number = int(matches_number)

    with open(sam_file, 'r') as sams:
        for sam in sams:
            if not sam.startswith('@'):
                full_read_name, sam_flag, chromosome, position, mapq, cigar, name_of_mate, position_of_mate, template_length, read_sequence, read_quality = sam.strip().split()[
                                                                                                                                                            :11]
                cigar_dic = split_string_to_dict(cigar)
                if 'M' in cigar_dic.keys():
                    sum_match = sum(cigar_dic['M'])
                else:
                    sum_match = 0
                if full_read_name not in dic.keys():
                    dic[full_read_name] = []
                    dic[full_read_name].append([position, sum_match])
                else:
                    dic[full_read_name].append([position, sum_match])
    for i in dic.keys():
        result = False
        max_match = 0
        for j in dic[i]:
            position = int(j[0])
            sum_match = int(j[1])
            if sum_match > matches_number:
                result = True
                if sum_match > max_match:
                    max_match = sum_match
            else:
                if (position + sum_match) > matches_number:
                    result = True
                if sum_match > max_match:
                    max_match = sum_match

        if result:
            outfile.write(i + '\t' + str(max_match) + '\n')

    outfile.close()
    
def split_string_S(s):
    numbers = re.findall(r'\d+', s)
    letters = re.findall(r'[a-zA-Z]+', s)

    if letters == ['S', 'M', 'S']:
        if int(numbers[0]) > int(numbers[2]):
            return 0, int(numbers[0])
        else:
            return int(numbers[0])  + int(numbers[1]), int(numbers[0])  + int(numbers[1]) + int(numbers[2])
    elif letters == ['S', 'M']:
        return 0, int(numbers[0])
    elif letters == ['M', 'S']:
        return int(numbers[0]), int(numbers[0])  + int(numbers[1])
    else:
        return False, False

def sam_split_S(sam_file, outfile, namefile):
    name_lis = []
    outfile = open(outfile, 'w')

    with open(namefile, 'r') as lls:
        for ll in lls:
            name_lis.append(ll.strip().split()[0])

    with open(sam_file, 'r') as sams:
        for sam in sams:
            if not sam.startswith('@'):
                full_read_name, sam_flag, chromosome, position, mapq, cigar, name_of_mate, position_of_mate, template_length, read_sequence, read_quality = sam.strip().split()[
                                                                                                                                                            :11]
                if full_read_name in name_lis:
                    if chromosome == '*':
                        seq_S = read_sequence
                        quality_S = read_quality
                        outfile.write('@' + full_read_name + '\n' + seq_S + '\n' + '+' + '\n' + quality_S + '\n')
                    else:
                        start, end = split_string_S(cigar)
                        if start or end:
                            seq_S = read_sequence[start:end]
                            quality_S = read_quality[start:end]
                            outfile.write('@' + full_read_name + '\n' + seq_S + '\n' + '+' + '\n' + quality_S + '\n')
    outfile.close()

def choose_max_matches(sam_file, outfile):
    outfile = open(outfile, 'w')
    dic = {}
    with open(sam_file, 'r') as sams:
        for sam in sams:
            if not sam.startswith('@'):
                full_read_name, sam_flag, chromosome, position, mapq, cigar, name_of_mate, position_of_mate, template_length, read_sequence, read_quality = sam.strip().split()[:11]
                cigar_dic = split_string_to_dict(cigar)
                if 'M' in cigar_dic.keys():
                    sum_match = sum(cigar_dic['M'])
                else:
                    sum_match = 0
                if full_read_name not in dic.keys():
                    dic[full_read_name] = [sum_match, sam.strip()]
                else:
                    if sum_match > dic[full_read_name][0]:
                        dic[full_read_name] = [sum_match, sam.strip()]
    for i in dic.keys():
        outfile.write(dic[i][1] + '\n')

    outfile.close()

def main():
    parameter = sys.argv[1]
    if parameter == 'guideseq_sam_positive':
        sam_file, outfile, matches_number = sys.argv[2:5]
        guideseq_sam_positive(sam_file, outfile, matches_number)

    elif parameter == 'sam_split_S':
        sam_file, outfile, namefile = sys.argv[2:5]
        sam_split_S(sam_file, outfile, namefile)
        
    elif parameter == 'choose_max_matches':
        sam_file, outfile = sys.argv[2:4]
        choose_max_matches(sam_file, outfile)

if __name__ == '__main__':
    main()
