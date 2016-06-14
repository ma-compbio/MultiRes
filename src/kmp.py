#!/usr/bin/env
import sys
def kmpTable(string):
    T = [-1]
    T.append(0)
    positionT = 2
    positionS = 0
    while positionT < len(string):
        if string[positionT-1] == string[positionS]:
            positionS += 1
            T.append(positionS)
            positionT += 1
        elif positionS > 0:
            positionS = T[positionS]
        else:
            T.append(0)
            positionT += 1
    return T

def kmp(string0,string1,T):
    match0 = 0
    match1 = 0
    match_points = []
    while match0+match1 < len(string0) and match1 < len(string1):
    #    print match0,match1,string0[match0],string1[match1]
        if string1[match1] == string0[match0+match1]:
            if match1 == len(string1)-1:
                match_points.append(match0)
                match1 = 0
                match0 += 1
            else:
#                match0 += 1
                match1 += 1
        else:
            if T[match1] > -1:
                match0 += match1-T[match1]
                match1 = T[match1]
            else:
                match1 = 0
                match0 += 1
    return match_points    
    
if __name__ == "__main__":
    string0 = sys.argv[1]
    string1 = sys.argv[2]

    print string0
    print string1
    T = kmpTable(string1)
    print kmp(string0,string1,T)
    
