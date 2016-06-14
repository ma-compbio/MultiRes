#!/usr/bin/env

import sys

def readProfiles(filename):
    cn_dict = {}
    with open(filename) as cn_file:
        for line in cn_file:
            mark,cn = line.strip().split(' ')
            cn_dict[mark] = int(cn)
    return cn_dict

def filterFamilies(filename,cn_dict,threshold):
    to_print = False
    with open(filename) as families_file:
        for line in families_file:
            if not(line.isspace()) and line[0] == '>':
                ident = line.strip()
                if cn_dict[ident] < threshold+1:
                    to_print = True
                else:
                    to_print = False
            if to_print:
                print line.strip()


if __name__ == "__main__" :
    copy_numbers = readProfiles(sys.argv[2])
    threshold = int(sys.argv[3])
    if threshold == -1:
        threshold = sys.maxint
    filterFamilies(sys.argv[1],copy_numbers,threshold)
