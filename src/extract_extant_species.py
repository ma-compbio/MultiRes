#!/usr/bin/env
import sys
import DataStructures.nhxtree as nhx

def treeClimber(T,v,marker = '@'):
    p = v.parent
    while p != None and p.node_name != marker:
        v = p
        p = v.parent
    return p

def printIngroupPairs(T,ingroups,marker,outfile):
    for i in range(len(ingroups)-1):
        for j in range(i+1,len(ingroups)):
            lca = nhx.queryLCA(T,[ingroups[i],ingroups[j]])
            if lca.node_name == marker:
                outfile.write(ingroups[i].node_name+' '+ingroups[j].node_name+'\n')

def printPairs(T,ing,outg,marker,to_filename):
    outfile = open(to_filename,'w')
    outfile.write("#outgroup\n")
    for i in ing:
        for o in outg:
            outfile.write(i.node_name+' '+o.node_name+"\n")
    outfile.write("#ingroup\n")
    printIngroupPairs(T,ing,marker,outfile)
    outfile.close()

def printExtantData(T,pairs_filename,marker='@'):
    leaves = T.leaves()
    ingroups = []
    outgroups = []
    for l in leaves:
        v = l
        if v.node_name != '':
            p = treeClimber(T,v,marker)
            if p != None and p.node_name == marker:
                ingroups.append(l)
            else:
                outgroups.append(l)
    printPairs(T,ingroups,outgroups,marker,pairs_filename)
    print "#INGROUPS"
    print '\n'.join([i.node_name for i in ingroups])
    print "#OUTGROUPS"
    print '\n'.join([o.node_name for o in outgroups])


if __name__ == "__main__":
    T = nhx.parse_nhx(open(sys.argv[1]).readline().strip())
    T.set_parent_root()
    printExtantData(T,sys.argv[2],'@')
