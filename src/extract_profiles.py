#!/usr/bin/env

import sys
import DataStructures.nhxtree as nhx

def getAncestralProfile(T,extant_profiles):
    max_prof = max([extant_profiles[s] for s in extant_profiles])
    ancestral_profiles = nhx.sankoff_rousseau(T,extant_profiles,range(max_prof+1))
    subtree = T.search('@')[0]
    return ancestral_profiles[subtree]


def printProfiles(filename,sp_names,tree):
    ident = ''
    with open(filename) as family_file:
        for line in family_file:
            if not(line.isspace()):
                if line[0] == '>':
                    if ident != '':
                        print ident + ' ' + str(getAncestralProfile(tree,new_profiles))
                        #print '\n'.join([x+' '+str(new_profiles[x]) for x in new_profiles])
                    new_profiles = dict([(x,0) for x in sp_names])
                    ident = line.strip()
                else:
                    current_sp = line.strip().split('.')[0]
                    new_profiles[current_sp] += 1
    if ident != '':
        print ident + ' ' + str(getAncestralProfile(tree,new_profiles))
        #print ident
        #print '\n'.join([x+' '+str(new_profiles[x]) for x in new_profiles])


if __name__ == "__main__":
    T = nhx.parse_nhx(open(sys.argv[1]).readline().strip())
    leaves = [x.node_name for x in T.leaves()]
    printProfiles(sys.argv[2],leaves,T)
