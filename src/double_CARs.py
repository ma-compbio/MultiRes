#!/usr/bin/env

import sys
from marker_string import *

def orientedEnds(marker):
    if marker[0] == '-':
        tail,head = readUndoubledMarker(marker[1:])
        return [head,tail]
    elif marker[0] == '+':
        tail,head = readUndoubledMarker(marker[1:])
        return [tail,head]
    else:
        tail,head = readUndoubledMarker(marker)
        return [tail,head]

def doubleCAR(content):
    new_content = []
    for c in content:
        new_content += orientedEnds(c)
    print "_Q "+' '.join(new_content)+" Q_"

def doubleCARs(filename):
    with open(filename) as carfile:
        for line in carfile:
            if not(line.isspace()):
                if line[0] == '>' or line[0] == '#':
                    print line.strip()
                else:
                    doubleCAR(line.strip().split(' ')[1:-1])


if __name__ == "__main__":
    doubleCARs(sys.argv[1])

