# encoding:utf-8

import sys

def file_to_lis(file):
    lis = []
    with open(file, 'r') as lls:
        for ll in lls:
            name = ll.strip()
            lis.append(name)
    return lis

def file_inter_file(file1, file2, not_in_file):
    lis1 = file_to_lis(file1)
    lis2 = file_to_lis(file2)
    not_in_file = open(not_in_file, 'w')
    result = list(set(lis1).difference(set(lis2)))
    for i in result:
        not_in_file.write(i + '\n')
    not_in_file.close()



def main():

    file1, file2, not_in_file = sys.argv[1:4]
    file_inter_file(file1, file2, not_in_file)

if __name__ == '__main__':
    main()
