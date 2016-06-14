#!/usr/bin/env

import sys

def convertCAR(line,counter):
    content = line.strip().split(' ')
    if '$' in content[-1]:
        print "#CAR"+str(counter)
        print "_Q " + ' '.join(content[:-1]) + " Q_\n"
    else:
        print "#CAR"+str(counter)
        print ' '.join(content) + '\n'
    counter += 1
    return counter

def convertCARs(filename):
    counter = 1
    with open(filename) as car_file:
        for line in car_file:
            if not(line.isspace()):
                if line[0] == '>':
                    continue
#                    print line.strip()
                elif line[0] != '#':
                    counter = convertCAR(line,counter)


if __name__ == "__main__" :
    convertCARs(sys.argv[1])
