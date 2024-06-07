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


def hifi_sam_classification(sam_file, negative_chrname, matches_number, outfile, unknown_file, positive_file):
    outfile = open(outfile, 'w')
    unknown_file = open(unknown_file, 'w')
    positive_file = open(positive_file, 'w')
    matches_number = int(matches_number)
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
                    dic[full_read_name] = {}
                    dic[full_read_name][chromosome] = sum_match
                else:
                    if chromosome not in dic[full_read_name].keys():
                        dic[full_read_name][chromosome] = sum_match
                    else:
                        max_sorce = max(sum_match, dic[full_read_name][chromosome])
                        dic[full_read_name][chromosome] = max_sorce
    for i in dic.keys():
        max_value = float('-inf')
        max_key = None
        result = False

        for inner_key, value in dic[i].items():
            if value > max_value:
                max_value = value
                max_key = inner_key
        # print(i, max_key, max_value)

        key_len = len(set(dic[i].keys()))

        if key_len == 1:
            result = True
        # print(key_len, dic[i].keys(), result)
        if result:
            if max_key == negative_chrname and max_value > matches_number:
                outfile.write(i + '\n')
            elif max_key == negative_chrname and max_value <= matches_number:
                unknown_file.write(i + '\n')
            else:
                positive_file.write(i + '\n')
        else:
            positive_file.write(i + '\n')

    outfile.close()
    unknown_file.close()
    positive_file.close()


def sam_extract(sam_file, name_file, out_file):
    out_file = open(out_file, 'w')
    name_lis = []

    with open(name_file, 'r') as lls:
        for ll in lls:
            name_lis.append(ll.strip())

    with open(sam_file, 'r') as sams:
        for sam in sams:
            if not sam.startswith('@'):
                full_read_name, sam_flag, chromosome, position, mapq, cigar, name_of_mate, position_of_mate, template_length, read_sequence, read_quality = sam.strip().split()[:11]
                if full_read_name in name_lis:
                    out_file.write(sam)

    out_file.close()

def main():

    parameter = sys.argv[1]
    if parameter == 'hifi_sam_classification':
        sam_file, negative_chrname, matches_number, outfile, unknown_file, positive_file = sys.argv[2:8]
        hifi_sam_classification(sam_file, negative_chrname, matches_number, outfile, unknown_file, positive_file)
    elif parameter == 'sam_extract':    
        sam_file, name_file, out_file = sys.argv[2:5]
        sam_extract(sam_file, name_file, out_file)


if __name__ == '__main__':
    main()