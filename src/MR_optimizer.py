#!/usr/bin/env

import sys
import copy
from DataStructures import markers
from DataStructures import genomes
from marker_string import *
from adjtools import *
import DataStructures.nhxtree as nhx
import networkx as nx
import linearization as lz
import kmp
import lcs
import weight as wt
import random

def getNeighbours(G,u,v,adjs):
    to_add = []
    for x in G[v]:
        if not(x in G[u]):
            if ((x[0],v[0]) in adjs or (v[0],x[0]) in adjs) and\
                    adjs[(x[0],v[0])] > G[v][x]['weight']:
                to_add.append((x,adjs[(v[0],x[0])]))
            else:
                to_add.append((x,G[v][x]['weight']))
    return to_add

# Contract the edge {u,v} to u, removing v and its mate.
def contractEdge(G,u,v,adjs):
    if u != v:
        to_add = getNeighbours(G,u,v,adjs)
        to_mate_add = getNeighbours(G,makeMate(u),makeMate(v),adjs)
        for a in to_add:
            x,weight = a
            G.add_edge(u,x,weight=weight)
        for a in to_mate_add:
            x,weight = a
            G.add_edge(makeMate(u),x,weight=weight)

def addToGeneLocusParent(gene_assoc,mark,loci_list):
    for l,c in loci_list:
        addToGeneLocus(gene_assoc,mark,c,l)


def addToGeneLocus(gene_assoc,mark,car_id,locus):
    if mark not in gene_assoc:
        gene_assoc[mark] = {}
    if car_id not in gene_assoc[mark]:
        gene_assoc[mark][car_id] = []
    if locus not in gene_assoc[mark][car_id]:
        gene_assoc[mark][car_id].append(locus)


# Multiplicities computed by Sankoff-Rousseau
def getMult(T,mults,printer=False):
    if printer:
        for m in mults:
            print "Infer function",m, mults[m]
    leaves = T.leaves()
    leaf_labels = {}
    possibilities = []
    for l in leaves:
        if l.node_name in mults:
            leaf_labels[l.node_name] = mults[l.node_name]
            if mults[l.node_name] not in possibilities:
                possibilities.append(mults[l.node_name])
        else:
            leaf_labels[l.node_name] = 0#mults[l.node_name]
            if 0 not in possibilities:
                possibilities.append(0)
    new_mult = nhx.sankoff_rousseau(T,leaf_labels,possibilities)
    return new_mult


def findTandems(adjs):
    tandems = {}
    for a in adjs:
        a0,a1 = a
        if a0 == makeMate(a1):
            tandems[a0] = 1
            tandems[a1] = 1
        elif a0 == a1:
            tandems[a0] = -1
            tandems[a1] = -1
    return tandems


def findEvenCycles(adj_dict,component):
    visited = set([])
    is_cyclic = {}
    build_traversal = []
    stack = component[0]
    visited.add(component[0])
    while len(adj_dict[stack[-1]]) == 0:
        rand_pos = random.choice(component)
        stack = [rand_pos]
        visited.add(rand_pos)
    while stack != []:
        current = stack.pop()
        build_traversal.append(current)
        nbhood = adj_dict[current]
        for v in nbhood:
            if makeMate(v) not in visited:
                stack.append(makeMate(v))
            else:
                is_cyclic[v] = current



def getPathNew(adjs,mults,weighted_tuples,orientation={}):
    paths = []
    ends = getEndsNew(adjs,mults,orientation)
    used = {}
    adj_dict = tupleToDict(adjs,False)
    tandem_dups = findTandems(adjs)
    for i in range(len(ends)):
        current_ends = ends[i][:]
        current_path_set = []
        while current_ends != []:
            current_path = [current_ends.pop()]
            prev = current_path[-1]
            current_path.append[makeMate(prev)]
            current = current_path[-1]
            mults[prev] -= 1
            mults[current] -= 1
            while current != '' and max([mults[x] for x in mults]) > 0:
                if current in tandem_dups and\
                        len(adj_dict[current])-1 < mults[current]:
                    dir = tandem_dups[current]
                    if dir > 0:
                        prev = makeMate(current)
                    elif dir < 0:
                        prev = current
                else:
                    if len(adj_dict[current]) == 1:
                        prev = adj_dict[current][-1]
                    elif len(adj_dict[current]) > 1 and current in tandem:
                        nbhood = adj_dict[current]
                    elif len(adj_dict[current]) > 1:
                        for e in adj_dict[current]:
                            if e != current and e != makeMate(current):
                                if e not in ends:
                                    ends.append(e)
                current_path.append(prev)
                current = makeMate(prev)
                current_path.append(current)
                mults[prev] -= 1
                mults[current] -= 1
            current_path_set.append(current_path)                
        paths.append(mergePaths(current_path_set,adj_dict))


def mergePaths(path_list,adj_dict):
    new_path_list = []
    merged_paths = []
    for i in range(len(path_list)):
        e0,e1 = path_list[i][0],path_list[i][-1]
        if e0 in adj_dict or e1 in adj_dict:
            if e0 in adj_dict and e1 in adj_dict:
                append_dict = extendPath(path_list[i],i,0,len(path_list[i]),path_list,adj_dict)
            elif e0 in adj_dict:
                append_dict = extendPath(path_list[i],i,0,len(path_list[i]),path_list,adj_dict)
            elif e1 in adj_dict:
                append_dict = extendPath(path_list[i],i,0,len(path_list[i]),path_list,adj_dict)
                

def extendPath(seed,index,e0=-1,e1=-2,path_list=[],adj_dict={}):
    new_path = seed
    to_append = {}
    if e0 >= 0:
        for i in range(len(path_list)):
            if i != index:
                f0,f1 = path_list[i][0],path_list[i][-1]
                if f0 in adj_dict[seed[e0]]:
                    to_append[i] = (-1,1)
                    break                    
                elif f1 in adj_dict[seed[e0]]:
                    to_append[i] = (1,1)
                    break                    
    if e1 > -2:
        for i in range(len(path_list)):
            if i != index:
                f0,f1 = path_list[i][0],path_list[i][-1]
                if f0 in adj_dict[seed[e1-1]]:
                    to_append[i] = (1,-1)
                    break                    
                elif f1 in adj_dict[seed[e1-1]]:
                    to_append[i] = (-1,-1)
                    break                    
    return to_append


# Computing a `good' traversal of the graph.

def getPath(G,adjs,vertices,segment,gene_assoc,\
        perm_mults,mults,weighted_tuples,orientation={}):
    paths = []
    sys.stdout.flush()        
    # Set of ends, defined as markers 
    # with no neighbours
    ends = getEndsNewer(G,adjs,mults,orientation,vertices,segment,gene_assoc)
    used = {}
    for a in adjs:
        used[a] = 0
    weights = {}
    for w in weighted_tuples:
        if (w[0],w[1]) in adjs:
            weights[(w[0],w[1])] = w[-1]
        elif (w[1],w[0]) in adjs:
            weights[(w[1],w[0])] = w[-1]            
    path = []
    tandem_used = {}
    # Sample from ends, till no ends left
    while ends != []:
        current = ends.pop()
        path = []
        while current != '' and max([mults[x] for x in mults]) > 0:
            if mults[current] <= 0 or mults[makeMate(current)] <= 0:
                break
            path.append(current)
            path.append(makeMate(current))
            if current in mults:
                mults[current] -= 1
            else:
                mults[current] = mults[makeMate(current)]-1
            mults[makeMate(current)] -= 1 
            candidates = []
            is_tandem = -1
            for a in adjs:
                if makeMate(current) in a:
                    candidates.append(a)
                    if current in a and a not in tandem_used:
                        is_tandem = len(candidates)-1
                        tandem_used[a] = False
            prev = current
            for c in candidates:
                if is_tandem >= 0 and \
                        tandem_used[candidates[is_tandem]] == False and \
                        c != candidates[is_tandem]:
                    continue
                elif is_tandem >= 0 and tandem_used[candidates[is_tandem]] == False \
                        and c == candidates[is_tandem]:
                    tandem_used[c] = True
                c0,c1 = c
                prospect = ''
                if c0 == makeMate(current):
                    prospect = c1
                else:
                    prospect = c0
                if mults[prospect] == 0:
                    continue
                current = prospect
                break
            if prev == current and is_tandem < 0:
                current = ''
        paths.append(path)
    to_add = []
    for m in vertices:
        if mults[m] == perm_mults[m]:
            to_add.append(m)
    used_marks = dict([(x,False) for x in to_add])
    for v in to_add:
        if not(used_marks[v]):
            paths.append([v,makeMate(v)])
            used_marks[v] = True
            used_marks[makeMate(v)] = True
            mults[v] -= 1
            mults[makeMate(v)] -= 1
    return paths        

# Calculate where to start the traversal
def getEnds(adjs,mults,orientations):
    # Ends list
    ends = []
    # Counts checks if number of adjacencies 
    # is strictly less than number of edges.
    # If so, possibility for end.
    counts = {}
    # Iterate over optimized adjacencies
    for i in range(len(adjs)-1):
        a0,a1 = adjs[i]
        # Assume we want to check if 
        # their mates are ends
        counts[makeMate(a0)] = 0
        counts[makeMate(a1)] = 0
        # Reiterate over all adjacencies
        for j in range(len(adjs)):
            b0,b1 = adjs[j]
#            if b0 == makeMate(b1):
#                continue
            if j == i:
                continue
            # The mate(s) are found in adjacencies
            if makeMate(a0) in adjs[j]:
                if makeMate(a0) not in counts:
                    counts[makeMate(a0)] = 0
                counts[makeMate(a0)] += 1
            if makeMate(a1) in adjs[j]:
                if makeMate(a1) not in counts:
                    counts[makeMate(a1)] = 0
                counts[makeMate(a1)] += 1
            # No spare copies left    
            if counts[makeMate(a0)] == mults[makeMate(a0)]\
                    and counts[makeMate(a1)] == mults[makeMate(a1)]:
                break
        # Extra! Extra! Read all about it!
        # In other words, there are spare copies to act as ends.
        if counts[makeMate(a0)] < mults[makeMate(a0)]:
            ends.append(makeMate(a0))
        if counts[makeMate(a1)] < mults[makeMate(a1)]:
            ends.append(makeMate(a1))
    if ends == []:
        connected_components = findConnectedComponents(adjs)
        for c in connected_components:
            ends.append(random.choice(c))
    return ends



def getEndsNewer(G,adjs,mults,orientations,\
        vertices,segment,gene_assoc):
    a_dict = {}
    for e in adjs:
        if e[0] not in a_dict:
            a_dict[e[0]] = []
        a_dict[e[0]].append(e[1])    
        if e[1] not in a_dict:
            a_dict[e[1]] = []
        a_dict[e[1]].append(e[0])    
    ends = []
    for v in G:
        if makeMate(v) in a_dict and \
                v not in a_dict:
            ends.append(v)                    
    if len(ends) == 0:
        print vertices        
#        getEndsInner(G,mults,segment,gene_assoc)
    return ends            
#    doStuff()

def getEndsNew(G,adjs,mults,orientations,\
        vertices,segment,gene_assoc):
    # Ends list
    ends = []
    # Counts checks if number of adjacencies 
    # is strictly less than number of edges.
    # If so, possibility for end.
    counts = {}
    connected_components= findConnectedComponents(adjs)
#    print "CC:",len(connected_components)
    adj_dict = tupleToDict(adjs,False)
    # Iterate over optimized adjacencies
    for c in connected_components:
        current_ends = []
        for a in c:
            if len(adj_dict[a]) == 0:
                current_ends.append(a)
        if len(current_ends) == 0:
            for a in c:
                nbhood = adj_dict[a]
                if len(adj_dict[a]) < mults[a] and\
                        len(adj_dict[makeMate(a)]) > len(adj_dict[a]):
                    current_ends.append(a)            
        if len(current_ends) == 0:
            current_ends.append(c[-1])
        ends.append(current_ends)
    ends = [x[0] for x in ends]        
    return ends


# Given an interval and the profiles, 

def findEdge(interval,profiles):
    for i in range(len(interval)):
        x,y = interval[i]
        if makeMate(x) == y and profiles[x] == 1:
            for u,v in interval:
                if (u,v) != (x,y) and (u,v) != (y,x) and ( u == x or v == x):
                    return i,0
                elif (u,v) != (x,y) and (u,v) != (y,x) and ( u == y or v == y):
                    return i,y
    return -1,None        

# Find the parsimonious content of the interval.

def findParsimoniousContent(locus,nhx_file,species,to_test,sp_mults,profiles):
    # locus is block in low resolution.
    # nhx_file is the species tree
    # species is the set of adjacencies between high res blocks in 
    #      each extant species
    # to_test is the total set of all adjacencies.
    # sp_mults is the number of local copies of a high res block in a given species
    # Also keep track of all vertices encountered. It is possible that they are spread out.
    # Call this V
    T = None
    with open(nhx_file) as tree_file:
        for line in tree_file:
            if not(line.isspace()):
                T = nhx.parse_nhx(line.strip())
    T.set_parent_root()
    leaves = T.leaves()
    keep = []
    species_dict = {}
    total_sp = []
    for x in to_test:
        current_labels = {}
        for s in species:
            if s not in total_sp:
                total_sp.append(s)
            for chr in species[s]:
                if x in species[s][chr] or tuple(x[::-1]) in species[s][chr]:
                    if x not in species_dict:
                        species_dict[x] = []
                    species_dict[x].append(s)
                    current_labels[s] = '1'
            if s not in current_labels:
                current_labels[s] = '0'
        new_labels = T.dollo_parsimony(current_labels,'1')
        need_tree = T.search('@')[0]
        actual_label = new_labels[need_tree]
        if actual_label == '1':
            keep.append(x)
    marker_mults = {}
    for x in sp_mults:
        printer = False
        need_tree = T.search('@')[0]
        new_mults = getMult(T,sp_mults[x],printer)
        marker_mults[x] = int(new_mults[need_tree])
        marker_mults[makeMate(x)] = int(new_mults[need_tree])
    zero_markers = [x for x in marker_mults if marker_mults[x] == 0]
    to_remove = []
    for i in range(len(keep)):
        a,b = keep[i]
        if marker_mults[a] == 0 or marker_mults[b] == 0:
            to_remove.append(i)
    for x in sorted(to_remove, reverse=True):
        keep.pop(x)
    edge_weights = []
    edge_weights = wt.weight_intervals_local(keep,species_dict,nhx_file)
    kept_tuples = [(keep[i][0],keep[i][1],edge_weights[i]) for i in range(len(edge_weights))]
    for x in zero_markers:
        del marker_mults[x]
    return kept_tuples,marker_mults,species_dict



def addToMarkerList(marker,markers,counter,labelling,local_mults,global_mults):
    if marker not in counter and marker not in labelling:
        markers[marker] = [(marker,0)]
        markers[makeMate(marker)] = [(makeMate(marker),0)]
        counter[marker] = 0
        counter[makeMate(marker)] = 0
        global_mults[(marker,0)] = local_mults[marker]
        global_mults[makeMate((marker,0))] = local_mults[marker]
        labelling[marker] = (marker,0)
        labelling[makeMate(marker)] = (makeMate(marker),0)
    elif marker not in labelling:
        counter[marker] += 1
        counter[makeMate(marker)] += 1
        markers[marker].append((marker,counter[marker]))
        markers[makeMate(marker)] = [(makeMate(marker),counter[makeMate(marker)])]
        global_mults[(marker,counter[marker])] = local_mults[marker]
        global_mults[makeMate((marker,counter[marker]))] = local_mults[marker]
        labelling[marker] = (marker,counter[marker])
        labelling[makeMate(marker)] = (makeMate(marker),counter[marker])
#    return newly_added

# Add new nodes and edges to the adjacency graph

def augmentGraph(G,markers,adjs,local_mults,global_mults,counter):
    added = []
    new_label = {}
    for e in adjs:
        u,v,weight = e
        if u not in new_label or makeMate(u) not in new_label:
            addToMarkerList(u,markers,counter,new_label,local_mults,global_mults)
            G.add_node(new_label[u])
            G.add_node(makeMate(new_label[u]))
        if v not in new_label or makeMate(v) not in new_label:            
            addToMarkerList(v,markers,counter,new_label,local_mults,global_mults)
            G.add_node(new_label[v])
            G.add_node(makeMate(new_label[v]))
        G.add_edge(new_label[v],new_label[u],weight=weight)
    for v in local_mults:
        if v not in new_label or makeMate(v) not in new_label:
            addToMarkerList(v,markers,counter,new_label,local_mults,global_mults)
            G.add_node(new_label[v])
            G.add_node(makeMate(new_label[v]))
    return new_label

# Add unaccounted for edges that are parsimoniously conserved

def addAdjs(G,c0,c1,adjs,marker_assoc,gene_assoc,global_mults,used,merge_dict={}):
    to_remove = []
    for v0 in marker_assoc[c0]:
        to_merge = False
        old = None
        for v1 in marker_assoc[c1]:
            if v0 == v1:
                if marker_assoc[c1][v1] != marker_assoc[c0][v0]:
                    to_merge = True
                    break
            elif (v0,v1) in adjs or (v1,v0) in adjs and (v0,v1) not in used and (v1,v0) not in used:
                G.add_edge(marker_assoc[c0][v0],marker_assoc[c1][v1],weight=1*adjs[(v0,v1)])
                used.add((v0,v1))
                used.add((v1,v0))
        if to_merge:
            merge_vertex = marker_assoc[c0][v0]
            if marker_assoc[c0][v0] in merge_dict:
                merge_vertex = merge_dict[marker_assoc[c0][v0]]
            contractEdge(G,merge_vertex,marker_assoc[c1][v1],adjs)
            original = marker_assoc[c1][v1]
            merge_dict[marker_assoc[c1][v1]] = merge_vertex 
            global_mults[merge_vertex] = max(global_mults[merge_vertex],\
                    global_mults[original])
            global_mults[makeMate(merge_vertex)] = max(global_mults[makeMate(merge_vertex)],\
                    global_mults[makeMate(original)])
            for car_id in gene_assoc[original]:
                for l in gene_assoc[original][car_id]:
                    marker_assoc[l].update({v0:merge_vertex})
                    marker_assoc[l].update({makeMate(v0):makeMate(merge_vertex)})
                    addToGeneLocus(gene_assoc,merge_vertex,car_id,l)
                    addToGeneLocus(gene_assoc,makeMate(merge_vertex),car_id,l)
#            for l in marker_assoc:
#                if v1 in marker_assoc[l]:
#                    gene_assoc[(v0,marker_assoc[c0][v0][-1])] += gene_assoc[(v1,marker_assoc[c1][v1][-1])]
            to_remove.append(original)
    return to_remove,used 
#            G.remove_node(marker_assoc[c1][v1])        
#            G.remove_node(makeMate(marker_assoc[c1][v1]))


def addToWindow(G,v,marker_assoc,gene_assoc,block_assocs):
    current_stack = [v]
    assignable = []
    prev_len = -1
    to_add_to = []
    traversed = []
    while len(current_stack) != 0 and assignable == []:
        peeked_vertex = current_stack.pop()
        to_add_to.append(peeked_vertex)
        if peeked_vertex not in gene_assoc and \
                assignable == []:
            for u in G[peeked_vertex]:
                if u not in current_stack and \
                        u not in to_add_to and \
                        u not in gene_assoc:
                    current_stack.append(u)
                elif u in gene_assoc and \
                        u not in traversed:
                    traversed.append(u)
                    for car_id in gene_assoc[u]:
                        for l in gene_assoc[u][car_id]:
                            if l not in assignable:
                                assignable.append((l,car_id))
            peeked_mate = makeMate(peeked_vertex)
            if peeked_mate not in current_stack and \
                    peeked_mate not in to_add_to:
                current_stack.append(peeked_mate)
    if assignable != []:            
        for to_assign in to_add_to:
            current_mate = makeMate(to_assign)
            ident,locality = to_assign
            addToGeneLocusParent(gene_assoc,to_assign,assignable)
            addToGeneLocusParent(gene_assoc,current_mate,assignable)
            for l,c in assignable:
                if ident not in marker_assoc[l]:
                    marker_assoc[l][ident] = to_assign 
                    marker_assoc[l][makeMate(ident)] = to_assign
                else:
                    for x in G[to_assign]:
                        G.add_edge(marker_assoc[l][ident],x,weight=G[to_assign][x]['weight'])
                    for y in G[current_mate]:
                        G.add_edge(makeMate(marker_assoc[l][ident]),y,weight=G[current_mate][y]['weight'])

def returnAssociatedVertices(G,association_subdict):
    vertices = []
    keep_to_remove = []
    for v in association_subdict:
        vv = association_subdict[v]
        if vv not in vertices and vv in G:
            vertices.append(vv)
        elif vv not in G and \
                vv not in keep_to_remove:
            keep_to_remove.append(v)
    return vertices,keep_to_remove

def returnInducingVertices(G,association_dict, \
        consecutive_markers,car_id,gene_assoc):
    vertices_to_induce = []
    for i in range(len(consecutive_markers)-1):
        c0 = consecutive_markers[i]
        c1 = consecutive_markers[i+1]
        locus = tuple(sorted([c0,c1]))
        if locus in association_dict:
            to_induce,to_remove = returnAssociatedVertices(G,association_dict[locus])
            vertices_to_induce += to_induce
            for v in to_remove:
                if v in association_dict[locus]:
                    del association_dict[locus][v]
    addAdjacentNonInducedVertices(G,association_dict,vertices_to_induce,car_id,gene_assoc)
    vertices_to_induce = list(set(vertices_to_induce))
    return vertices_to_induce

def addAdjacentNonInducedVertices(G,association_dict,vertices_to_induce,car_id,gene_assoc):
    more_to_add = []
    for v in vertices_to_induce:
        #if v not in G:
        #    G.add_node(v)
        #    G.add_node(makeMate(v))
        #    continue
        for u in G[v]:
            if v not in gene_assoc:
                loci_list = [(l,car_id) for l in association_dict\
                        if v in [association_dict[l][x] for x in association_dict[l]]]
                addToGeneLocusParent(gene_assoc,v,loci_list)
                addToGeneLocusParent(gene_assoc,makeMate(v),loci_list)
            if car_id not in gene_assoc[v]:                
                loci_list = [(l,car_id) for l in association_dict\
                        if v in [association_dict[l][x] for x in association_dict[l]]]
                addToGeneLocusParent(gene_assoc,v,loci_list)
                addToGeneLocusParent(gene_assoc,makeMate(v),loci_list)
            if not(isInAssociation(u,association_dict,gene_assoc[v][car_id])):
                more_to_add.append(u)
                more_to_add.append(makeMate(u))
    vertices_to_induce += more_to_add

def isInAssociation(vertex,association_dict,list_of_loci):
    for l in list_of_loci:
        if vertex in association_dict[l]:
            return True
    return False


def constructSubgraph(G,association_dict,consecutive_markers,car_id,gene_assoc):
    vertex_subset = returnInducingVertices(G,association_dict,\
            consecutive_markers,car_id,gene_assoc)
    Gprime = G.subgraph(vertex_subset)
    return Gprime

def pathConstructor(G,opt_edges,segment,gene_assoc,copy_numbers,orientation):
    profiles = copy.deepcopy(copy_numbers)
    kept_weights = []
    for e in opt_edges:
        v0,v1 = e
        weight = G[v0][v1]['weight']
        kept_weights.append((v0,v1,weight))
    return getPath(G,opt_edges,G.nodes(),segment,gene_assoc,copy_numbers,profiles,kept_weights,orientation)

def naivePathConstructor(Gprime,kept_adjs):
    paths = []
    for e in kept_adjs:
        paths.append([makeMate(e[0]),\
                e[0],e[1],\
                makeMate(e[1])])
    return paths

def findSubCAR(G,association_dict,gene_assoc,car_id,\
        consecutive_markers,copy_numbers,orientation):
    Gprime = constructSubgraph(G,association_dict,consecutive_markers,car_id,gene_assoc)
    kept_adjs = lz.linearization(Gprime,copy_numbers,1)
    paths = naivePathConstructor(Gprime,kept_adjs)
    return paths

def findSubCARs(G,association_dict,gene_assoc,\
        big_cars,copy_numbers,segment_size,orientations={}):
    tiling = {}
    for c in big_cars:
        tiling[c] = tileCAR(big_cars[c],G,association_dict,c,gene_assoc,\
                copy_numbers, segment_size,orientations)
    return tiling

def tileCAR(car, G, association_dict, car_id, gene_assoc,\
        copy_numbers, segment_size,orientation):
    tiling_paths = []
    start = 0
    while start < len(car) and (start+segment_size < len(car)+1 or segment_size > len(car)):
        j = 1
        while start+j < len(car) and j < segment_size:
            j += 1
        current_content = car[start:start+j]
        tiling_paths.append(findSubCAR(G,association_dict,gene_assoc,car_id,\
                current_content,copy_numbers,orientation))
        start += 1
    return tiling_paths

def unlocalize(car):
    new_car = [x[0] if x.__class__ == (1,3).__class__ else x for x in car]
    return new_car

def balanceCopies(G,u,u_content,copy_numbers,profile,\
        loci_dict,gene_dict):
    u_m = makeMate(u)
    print_tag = False
    if u == '14' or u == '13':
        print_tag = False 
    deletable = []        
    for x in u_content:
        if (u,x) in G and \
                (len(G[(u,x)]) == 0 and len(G[(u_m,x)]) == 0) and\
                x not in deletable:
            current_maxer = x
            deletable.append(x)
        elif copy_numbers[(u,x)] == 0:
            continue
    for x in deletable:
        copy_numbers[(u,x)] = 0
        copy_numbers[(u_m,x)] = 0
    while profile != sum([copy_numbers[(u,k)] for k in u_content]):
        if sum([copy_numbers[(u,k)] for k in u_content]) > profile:
            current_maxer = None
            current_max = -1
            for x in u_content: 
                if copy_numbers[(u,x)] == 0:
                    continue
                elif copy_numbers[(u,x)] > current_max:
                    current_max = copy_numbers[(u,x)]
                    current_maxer = x
            copy_numbers[(u,current_maxer)] -= 1
            copy_numbers[(u_m,current_maxer)] -= 1
        else:            
            current_minner = None
            current_min_nbs = 0
            current_max_nbs = -1
            current_max_nbrr = None
            for x in u_content: 
                #if copy_numbers[(u,x)] == len(G[(u,x)]):
                #    continue
                if copy_numbers[(u,x)]-len(G[(u,x)]) < current_min_nbs :
                    current_min_nbs = copy_numbers[(u,x)]-len(G[(u,x)])
                    current_minner = x
                if current_max_nbs < len(G[(u,x)]):
                    current_max_nbrr = x
                    current_max_nbs = len(G[(u,x)])
            if current_minner != None:        
                copy_numbers[(u,current_minner)] += 1
                copy_numbers[(u_m,current_minner)] += 1
            else:
                copy_numbers[(u,current_max_nbrr)] += 1
                copy_numbers[(u_m,current_max_nbrr)] += 1
    for x in u_content:
        if copy_numbers[(u,x)] == 0:
            G.remove_node((u,x))
            G.remove_node(makeMate((u,x)))
            for car_id in gene_dict[(u,x)]:
                for l in gene_dict[(u,x)][car_id]:
                    if u in loci_dict[l]:
                        del loci_dict[l][u]
                    if u_m in loci_dict[l]:
                        del loci_dict[l][u_m]
            del gene_dict[(u,x)]
            del gene_dict[(u_m,x)]
            del copy_numbers[(u,x)]
            del copy_numbers[(u_m,x)]
    if print_tag:
        print u,profile
        for i in u_content:
            print (u,i), copynumbers[(u,i)]

def addExtraContent(G,car_id,car_seq,blocks,sp_contents,nhx_file,marker_assoc,gene_assoc,\
        small_window,segment_length,copy_numbers,block_orients,counter,pars_adjs):
    for i in range(0,len(car_seq)-small_window):
        j = i+1
        while j < len(car_seq) and j-i < small_window:
            j += 1
        current_window = car_seq[i:j]
        #extended_window = car_seq[max(i-(segment_length>>2),0):min(i+(segment_length>>2),len(car_seq))]
        subAddExtraContent(G,car_id,current_window,current_window,blocks,sp_contents,nhx_file,marker_assoc,gene_assoc,\
                copy_numbers,block_orients,counter,pars_adjs)

def findInSpecies(block_dict,l_index,locus,window_seq,previous_content,\
        block_orients,prev_orient,candidate_loci,anc_orient):
    block = block_dict[locus]
    sp_found = []
    # For each locus in the block, do
    for sp in block:
        # Mark that locus is present in species sp
        sp_found.append(sp)
        # If no path has ever been started, add to path dict
        if sp not in window_seq:
            window_seq[sp] = {}
            previous_content[sp] = {}
            prev_orient[sp] = {}
            #current_block_sp.append(sp) # species that have current_block
        for chrom in block[sp]:
            current_content = block[sp][chrom]
#            print current_content
            checkInChrom(block,l_index,locus,sp,chrom,window_seq,\
                    current_content,previous_content,\
                    block_orients,prev_orient,candidate_loci,anc_orient)
    # For each species that does not have the locus,
    # start a new path, since previous one cannot be extended
    for sp in window_seq:
        if sp not in sp_found:
            # It is not found in any chromosome
            for chrom in window_seq[sp]:
                window_seq[sp][chrom].append([])

# For each locus, add it as an expected locus to the set of 
# genes it contains
def addToCandidateLoci(locus_index,content,candidate_loci):
    for mark in content:
        if mark not in candidate_loci:
            candidate_loci[mark] = []
        if makeMate(mark) not in candidate_loci:
            candidate_loci[makeMate(mark)] = []
        candidate_loci[mark].append(locus_index)
        candidate_loci[makeMate(mark)].append(locus_index)

def checkInChrom(block,l_index,locus,sp,chrom,window_seq,\
        current_content,prev_sp_content,\
        block_orients,prev_orient,candidate_loci,anc_orient):
    if chrom not in window_seq[sp]:
        window_seq[sp][chrom] = [[]]
    if chrom not in prev_sp_content[sp]:        
        prev_sp_content[sp][chrom] = None
    if chrom not in prev_orient[sp]:        
        prev_orient[sp][chrom] = 0
    addToCandidateLoci(l_index,current_content,candidate_loci)
    extendPath(locus,sp,chrom,window_seq,\
            current_content,prev_sp_content,\
            block_orients,prev_orient,anc_orient)

def setPrevContent(prev_sp_content,prev_orient,locus,sp,chrom,anc_orient):
    if anc_orient < 0:
        prev_sp_content[sp][chrom] = locus[0]
    elif anc_orient > 0:
        prev_sp_content[sp][chrom] = locus[-1]

# Tandem repeat adder. Only adds opposite tandems, of type
# (A -A).
def conditionalAdd(content,window_seq,prev_sp_content,\
        prev_orient,block_orients,locus,sp,chrom,anc_orient):
    current_orient = block_orients[locus][sp][chrom]
    if anc_orient < 0:
        if current_orient > 0:
            window_seq[sp][chrom][-1] += content[::-1]# + window_seq[sp][chrom][-1]
        else:
            window_seq[sp][chrom][-1] += content # + window_seq[sp][chrom][-1]
    elif anc_orient > 0:            
        if current_orient > 0:
            window_seq[sp][chrom][-1] += content# + window_seq[sp][chrom][-1]
        else:
            window_seq[sp][chrom][-1] += content[::-1] # + window_seq[sp][chrom][-1]
        #if len(window_seq[sp][chrom][-1]) > 0 and \
        #        content != window_seq[sp][chrom][-1][:2]:
        #    window_seq[sp][chrom][-1] = window_seq[sp][chrom][-1] + content
    prev_orient[sp][chrom] = current_orient
    setPrevContent(prev_sp_content,prev_orient,locus,sp,chrom,anc_orient)

# Extend path using the current content. 
# First, check if a tandem repeat exists.
# If yes, and the tandem is opposite, i.e. (A -A), add.
# Else, i.e. (A A), err on the side of caution, and do not add.
def extendPath(locus,sp,chrom,window_seq,\
        content,prev_sp_content,\
        block_orients,prev_orient,anc_orient):
    if prev_orient[sp][chrom] == 0:
#        print "Here1"
        window_seq[sp][chrom][-1] += content
        prev_orient[sp][chrom] = block_orients[locus][sp][chrom]
        setPrevContent(prev_sp_content,prev_orient,locus,sp,chrom,anc_orient)
    else:
        if prev_sp_content[sp][chrom] in locus or makeMate(prev_sp_content[sp][chrom]) in locus:
            if len(content) == 2:
#                print "Here2"
                conditionalAdd(content,window_seq,prev_sp_content,\
                        prev_orient,block_orients,locus,sp,chrom,anc_orient)
            else:
                if anc_orient < 0:
#                    print "Here3"
                    if block_orients[locus][sp][chrom] > 0:
                        window_seq[sp][chrom][-1] += content[::-1] #+ window_seq[sp][chrom][-1]
                    else:
                        window_seq[sp][chrom][-1] += content #+ window_seq[sp][chrom][-1]
                elif anc_orient > 0: 
#                    print "Here4"
                    if block_orients[locus][sp][chrom] > 0:
                        window_seq[sp][chrom][-1] += content
                    else:
                        window_seq[sp][chrom][-1] += content[::-1]
        else:
            if window_seq[sp][chrom][-1] == []:
#                print "Here6"
                window_seq[sp][chrom][-1] += content
            else:
#                print "Stuff",content
                window_seq[sp][chrom].append([])
                window_seq[sp][chrom][-1] += content
        prev_orient[sp][chrom] = block_orients[locus][sp][chrom]
        setPrevContent(prev_sp_content,prev_orient,locus,sp,chrom,anc_orient)

def subAddExtraContent(G,car_id,window,segment,blocks,parsed_adjs,nhx_file,marker_assoc,gene_assoc,\
        copy_nums,block_orients,counter,pars_adjs):
    in_species = {}
    adj_assignments = {}
    local_mults = {} # local multiplicities in a window
    window_sp_seq = {} # dict of species: dict of chrom: list of disjoint paths; each path=list of consecutive loci
    prev_orient = {} # dict of spec: dict of chrom?
    prev_sp = {}
    # time_keeper = []
    prev_block_sp = []
    candidate_loci = {}
    for j in range(len(window)-1):
        win = tuple(sorted([window[j],window[j+1]]))
        if win[0] != window[j]:
            anc_orient = -1
        else:
            anc_orient = 1
        #current_block_sp = []
        if win in blocks:
            findInSpecies(blocks,j,win,window_sp_seq,\
                    prev_sp,block_orients,prev_orient,candidate_loci,anc_orient)
    for sp in window_sp_seq:
        for chrom in window_sp_seq[sp]:
#            print sp,chrom,window_sp_seq[sp][chrom]
            for content in window_sp_seq[sp][chrom]:
                if content == []:
                    continue
                #for content in contents:
                for i in range(len(content)-1):
                    if content[i] not in local_mults:
                        local_mults[makeMate(content[i])] = {}
                        local_mults[content[i]] = {}
                    if sp not in local_mults[content[i]]:
                        local_mults[makeMate(content[i])][sp] = 0
                        local_mults[content[i]][sp] = 0
                    #local_mults[makeMate(content[i])][sp] += 1
                    local_mults[content[i]][sp] += 1
                    if i % 2 == 1:
                        a = tuple(sorted((content[i],content[i+1])))
                        if '-1' not in a:
#                            print a
                            if a not in adj_assignments: 
                                adj_assignments[a] = []
                            if sp not in adj_assignments[a]:    
                                adj_assignments[a].append(sp)    
                if content[i+1] not in local_mults:
                    local_mults[makeMate(content[i+1])] = {}
                    local_mults[content[i+1]] = {}
                if sp not in local_mults[content[i+1]]:
                    local_mults[makeMate(content[i+1])][sp] = 0
                    local_mults[content[i+1]][sp] = 0
                local_mults[makeMate(content[i+1])][sp] += 1
                local_mults[content[i+1]][sp] += 1
    parsAddEdges(G,car_id,adj_assignments,nhx_file,marker_assoc,gene_assoc,\
            window,segment,copy_nums,local_mults,counter,candidate_loci,pars_adjs)

def parsAddEdges(G,car_id,candidate_adjacencies,nhx_file,marker_assoc,gene_assoc,\
        window,segment,copy_numbers,profiles,counter,new_loci_assoc,pars_adjs):
    T = None
    # time_keeper = []
    with open(nhx_file) as tree_file:
        for line in tree_file:
            if not(line.isspace()):
                T = nhx.parse_nhx(line.strip())
    T.set_parent_root()
    leaves = T.leaves()
    species_dict = copy.deepcopy(candidate_adjacencies)
    keep = []
    candidate_markers = {}
    candidate_windows = {}
    for a in candidate_adjacencies:
        if a not in pars_adjs:
            continue
        leaf_labels = dict([(x,'1') for x in candidate_adjacencies[a]])
        new_labels = T.dollo_parsimony(leaf_labels,'1')            
        need_tree = T.search('@')[0]
        actual_label = new_labels[need_tree]
        if (actual_label == '1'):
            keep.append(a)
        elif (actual_label == '0'):
            del species_dict[a]    
    edge_weights = []
    edge_weights = wt.weight_intervals_local(keep,species_dict,nhx_file)
    kept_tuples = [(keep[i][0],keep[i][1],edge_weights[i]) for i in range(len(edge_weights))]
    #print "Pars time0:", # time_keeper[-1]-# time_keeper[-2]
    already_added = []
    for a in kept_tuples:
        a0,a1,weight = a
        #print a0,a1,makeMate(a0),makeMate(a1)
        already_added += [a0,a1,makeMate(a0),makeMate(a1)]
        addMissedEdge(G,car_id,a0,a1,marker_assoc,gene_assoc,\
                window,copy_numbers,profiles,weight,T,counter,new_loci_assoc)
    for m in profiles:
        if m not in candidate_markers:
            candidate_markers[m] = [sp for sp in profiles[m] if profiles[m][sp] > 0]
        if makeMate(m) not in candidate_markers:
            candidate_markers[makeMate(m)] = [sp for sp in profiles[makeMate(m)] if profiles[makeMate(m)][sp] > 0]
    for i in range(len(segment)-1):
        l = tuple(sorted([segment[i],segment[i+1]]))
        if l in marker_assoc:
            for x in marker_assoc[l]:
                already_added.append(x)
                already_added.append(makeMate(x))
    for m in candidate_markers:
        if m not in already_added and makeMate(m) not in already_added:
            leaf_labels = dict([(x,'1') for x in candidate_markers[m]])
            new_labels = T.dollo_parsimony(leaf_labels,'0')
            need_tree = T.search('@')[0]
            actual_label = new_labels[need_tree]
#            if m == '80':
#                print "Here",actual_label
            if (actual_label > 0):
                addMissedVertex(G,car_id,m,marker_assoc,gene_assoc,\
                        window,copy_numbers,profiles,T,counter,new_loci_assoc)
            already_added.append(m)
            already_added.append(makeMate(m))
    #print "Pars time1:", # time_keeper[-1]-# time_keeper[-2]


def addMissedVertex(G,car_id,m,marker_assoc,gene_assoc,\
        window,copy_numbers,profiles,tree,counter,loci_assoc):
    #if m == '38' or makeMate(m) == '38':
    #    print "HereVertex",car_id
    if m in counter:
        #if m == '38' or makeMate(m) == '38':
        #    print "Found"
        counter[m] += 1
        counter[makeMate(m)] += 1
    else:
        counter[m] = 0
        counter[makeMate(m)] = 0
    G.add_node((m,counter[m]))
    G.add_node((makeMate(m),counter[makeMate(m)]))    
    #if m == '38' or makeMate(m) == '38':
    #    print "HereVertex",car_id,m,counter[m]
    t_profile = getMult(tree,profiles[m])
    need_tree = tree.search('@')[0]
    profile = int(t_profile[need_tree])
    copy_numbers[(m,counter[m])] = profile
    copy_numbers[makeMate((m,counter[m]))] = profile
    min_win = min(loci_assoc[m])
    max_win = max(loci_assoc[m])+1
    for i in range(min_win,max_win):
        locus = tuple(sorted([window[i],window[i+1]]))
        if locus not in marker_assoc:
            marker_assoc[locus] = {}
        if m not in marker_assoc[locus]:
            marker_assoc[locus][m] = (m,counter[m])
            marker_assoc[locus][makeMate(m)] = makeMate((m,counter[m]))
        addToGeneLocus(gene_assoc,(m,counter[m]),car_id,locus)
        addToGeneLocus(gene_assoc,makeMate((m,counter[m])),car_id,locus)



def addMissedEdge(G,car_id,m0,m1,marker_assoc,gene_assoc,\
        window,copy_numbers,profiles,weight,tree,counter,loci_assoc):
    assoc_markers0 = []
    assoc_markers1 = []
    skippable = False    
    for i in range(len(window)-1):
        locus = tuple(sorted([window[i],window[i+1]]))
        if locus in marker_assoc:
            for x in marker_assoc[locus]:
                if m0 == m1:
                    if marker_assoc[locus][x][0] == m0:
                        assoc_markers0.append(marker_assoc[locus][x])
                        assoc_markers1.append(marker_assoc[locus][x])
                elif marker_assoc[locus][x][0] == m0:# or x[0] == m1:
                    assoc_markers0.append(marker_assoc[locus][x])
                    for v_other in G[marker_assoc[locus][x]]:
                        if v_other[0] == m1:
                            skippable = True
                elif marker_assoc[locus][x][0] == m1:# or x[0] == m1:
                    assoc_markers1.append(marker_assoc[locus][x])
                    for v_other in G[marker_assoc[locus][x]]:
                        if v_other[0] == m0:
                            skippable = True
    if not(skippable) and \
            (len(assoc_markers0) == 0 or len(assoc_markers1) == 0):
        if len(assoc_markers0) == 0:
            if m0 not in counter:
                counter[m0] = 0
                counter[makeMate(m0)] = 0
            else:
                counter[m0] += 1
                counter[makeMate(m0)] +=1
            max0 = counter[m0]         
            G.add_node((m0,max0))
            G.add_node(makeMate((m0,max0)))
            t_profile = getMult(tree,profiles[m0])
            need_tree = tree.search('@')[0]
            profile = int(t_profile[need_tree])
            copy_numbers[(m0,max0)] = profile
            copy_numbers[makeMate((m0,max0))] = profile
            assoc_markers0.append((m0,max0))
            min_win = min(loci_assoc[m0])
            max_win = max(loci_assoc[m0])+1
            for i in range(min_win,max_win):
                locus = tuple(sorted([window[i],window[i+1]]))
                if locus not in marker_assoc:
                    marker_assoc[locus] = {}
                if m0 not in marker_assoc[locus]:
                    marker_assoc[locus][m0] = (m0,max0)
                    marker_assoc[locus][makeMate(m0)] = makeMate((m0,max0))
                addToGeneLocus(gene_assoc,(m0,max0),car_id,locus)                    
                addToGeneLocus(gene_assoc,makeMate((m0,max0)),car_id,locus)                    
        if len(assoc_markers1) == 0:        
            if m1 not in counter:
                counter[m1] = 0
                counter[makeMate(m1)] = 0
            else:
                counter[m1] += 1
                counter[makeMate(m1)] +=1
            max1 = counter[m1]         
            G.add_node((m1,max1))
            G.add_node(makeMate((m1,max1)))
            t_profile = getMult(tree,profiles[m1])
            need_tree = tree.search('@')[0]
            profile = int(t_profile[need_tree])
            copy_numbers[(m1,max1)] = profile
            copy_numbers[makeMate((m1,max1))] = profile
            assoc_markers1.append((m1,max1))
            min_win = min(loci_assoc[m1])
            max_win = max(loci_assoc[m1])+1
            for i in range(min_win,max_win):
                locus = tuple(sorted([window[i],window[i+1]]))
                if locus not in marker_assoc:
                    marker_assoc[locus] = {}
                if m1 not in marker_assoc[locus]:
                    marker_assoc[locus][m1] = (m1,max1)
                    marker_assoc[locus][makeMate(m1)] = makeMate((m1,max1))
                addToGeneLocus(gene_assoc,(m1,max1),car_id,locus)
                addToGeneLocus(gene_assoc,makeMate((m1,max1)),car_id,locus)                    
    #print assoc_markers0
    #print assoc_markers1
#    print "Marker add:",# time_keeper[-1]-# time_keeper[-2]
    if not skippable:
        need_tree = tree.search('@')[0]
        t_profile0 = getMult(tree,profiles[m0])
        profile0 = int(t_profile0[need_tree])
        t_profile1 = getMult(tree,profiles[m1])
        profile1 = int(t_profile1[need_tree])
        for (m0,max0) in assoc_markers0:         
            for (m1,max1) in assoc_markers1:
#            print (m0,max0), (m1,max1)
                if profile0 > copy_numbers[(m0,max0)]:
                    copy_numbers[(m0,max0)] = profile0
                    copy_numbers[makeMate((m0,max0))] = profile0
                if profile1 > copy_numbers[(m1,max1)]:
                    copy_numbers[(m1,max1)] = profile1
                    copy_numbers[makeMate((m1,max1))] = profile1
                G.add_edge((m0,max0),(m1,max1))
                G[(m0,max0)][(m1,max1)]['weight'] = weight
#    print "Mults and edge add:",# time_keeper[-1]-# time_keeper[-2]
#            if m1 == '364':
#                print "NBhood",G[(m1,max1)]

def findDist(gene_assoc,car_id,loc_list,nbhood,index,blk_indices):
    if car_id not in gene_assoc[nbhood[index]]:
        return sys.maxint
    new_loc_list = gene_assoc[nbhood[index]][car_id]
    new_index_list = [blk_indices[l][car_id] for l in new_loc_list]
    differ = 0
    for i in new_index_list:
        if max([abs(l-lprime) for lprime in new_index_list for l in loc_list]) > differ:
            far = index
            differ = max([abs(l-lprime) for lprime in new_index_list for l in loc_list])
    return differ

def deleteRedundantEdges(G,v,nbhood,gene_assoc,blk_indices):
    i = 0
    car_assoc = gene_assoc[v]
    loc_list = []
    deletable = []
    for c in car_assoc:
        loc_list = [blk_indices[l][c] for l in gene_assoc[v][c]]
        while i < len(nbhood):
            farthest = -1
            diff = 0
            j = i
            farthest_dict = {}
            while j < len(nbhood) and \
                    nbhood[j][0] == nbhood[i][0]:
                farthest_dict[nbhood[j]] =  findDist(gene_assoc,c,loc_list,nbhood,j,blk_indices)
                j += 1                
            deletable_candidates = sorted([nbhood[k] for k in range(i,j)],key=lambda x: farthest_dict[x])
            deletable += deletable_candidates[1:]
            i = j            
        for x in deletable:
            if x in G[v]:
                G.remove_edge(v,x)



def removeSubCARS(car_list,copies):
    to_remove = []
    for cc in car_list.keys():
        if cc not in car_list or cc in to_remove:
            continue
        orig_car = car_list[cc]
        temp_car = unlocalize(car_list[cc])
        removed = False
        for c in car_list.keys():
            if c not in car_list or c == cc or c in to_remove:
                continue
            localized_candidate = car_list[c]
            candidate = unlocalize(car_list[c])
            if len(candidate) > len(temp_car):
                removable = True
                uniq_count = 0
                for i in range(len(temp_car)):
                    if copies[orig_car[i]] == 1 and temp_car[i] in candidate:
                        uniq_count += 1
                        continue
                    elif copies[orig_car[i]] == 1 and temp_car[i] not in candidate:
                        removable = False
                        break
                if removable and uniq_count > len(candidate)>>1:
#                    print "AA",' '.join(candidate)
#                    print "AB",' '.join(temp_car)
                    removed = True
                    to_remove.append(cc)
                    break
                match_pos = kmp.kmp(candidate,temp_car)
                match_neg = kmp.kmp(candidate,temp_car[::-1])
                if len(match_pos) > 0 or len(match_neg) > 0 or\
                    lcs.longestSubstring(candidate,temp_car)[0] >= len(temp_car)>>1:
                    #print match_pos
                    #print match_neg
                    to_remove.append(cc)
                    removed = True
                    break
                elif lcs.localEditDistance(candidate,temp_car)[-1] < len(temp_car)/10 or\
                        lcs.localEditDistance(candidate,temp_car[::-1]) < len(temp_car)/10:
                    to_remove.append(cc)
                    removed = True
                    break
            else:
                if len(candidate) < temp_car:
                    removable = True
                    uniq_count = 0
                    for i in range(len(candidate)):
                        if copies[localized_candidate[i]] == 1 and candidate[i] in temp_car:
                            uniq_count += 1
                            continue
                        elif copies[orig_car[i]] == 1 and temp_car[i] not in candidate:
                            removable = False
                            break
                    if removable and uniq_count > len(candidate)>>1:
                        #print "BA",' '.join(candidate)
                        #print "BB",' '.join(temp_car)
                        to_remove.append(c)
                        continue
                opp_match_pos = kmp.kmp(temp_car,candidate)
                opp_match_neg = kmp.kmp(temp_car[::-1],candidate)
                if len(opp_match_pos) > 0 or len(opp_match_neg) > 0 or\
                    lcs.longestSubstring(candidate,temp_car)[0] >= len(candidate)>>1:
                    #print opp_match_pos
                    #print opp_match_neg
                    to_remove.append(c)
                elif lcs.localEditDistance(temp_car,candidate)[-1] < len(candidate)/10 or\
                        lcs.localEditDistance(temp_car,candidate[::-1]) < len(candidate)/10:
                    to_remove.append(c)
        if removed:
            continue
    if len(to_remove) > 0:
        for c in to_remove:
            del car_list[c]

tree_file = sys.argv[1]

car_dict = {}
num = 0
#car_adj_set = set([])
block_indices = {}
with open(sys.argv[4]) as car_file:
    for line in car_file:
        if not(line.isspace()) and \
                not(line[0] == '#') and \
                not(line[0] == '>'):
            car = line.strip().split(' ')[1:-1]                    
            car_dict[car_id] = car
            for i in range(len(car_dict[car_id])-1):
#                car_adj_set.add(tuple(sorted([car_dict[car_id][i],car_dict[car_id][i+1]])))
                block_indices[tuple(sorted([car_dict[car_id][i],car_dict[car_id][i+1]]))] = {}
                block_indices[tuple(sorted([car_dict[car_id][i],car_dict[car_id][i+1]]))][car_id] = i
            num += 1
        elif not(line.isspace()) and \
                line[0] == '#':
            car_id = line[line.find('CAR'):].strip()



locus = ''
in_species = {}
make_locus = {}
local_copies = {}
copynumbers = {}
kept_species = []
possible_content = []
possible_markers = []
marker_sp_dict = {}
chaining = {}
content = {}
unparsed_content = {} # Because adjacencies and vertices may be spread out
addable = False
adj_orients = {}
big_adj_orients = {}
with open(sys.argv[2]) as containment_file:
    to_write = False
    is_first = True
    to_do = False
    for line in containment_file:
        if not(line.isspace()) and '>' in line:
            if is_first:
                is_first = False
            else:
                if len(in_species) > 0 and addable:
                    if locus not in content:
                        content[locus] = findParsimoniousContent(locus,tree_file,in_species,possible_content,local_copies,copynumbers)
                        unparsed_content[locus] = list_in_species 
                    else:
                        temp_content = findParsimoniousContent(locus,tree_file,in_species,possible_content,local_copies,copynumbers)
                        unparsed_content[locus] += list_in_species 
                        content[locus] += temp_content[0]
                        for k in temp_content[1]:
                            if k in content[locus][1] and \
                                    content[locus][1][k] < temp_content[1][k]:
                                content[locus][1][k] = temp_content[1][k]
                            elif k not in content[locus][1]:
                                content[locus][1][k] = temp_content[1][k]
                        for m in temp_content[2]:
                            if m not in content[locus][2]:
                                content[locus][2] = temp_content[2][m]
                            else:
                                for s in temp_content[2][m]:
                                    if s not in content[locus][2][m]:
                                        content[locus][2][m].append(m)
            locus = tuple(sorted(line[1:].strip().split(' ')))
            big_adj_orients[locus] = {}
            if len(locus) == 1:
                locus = tuple(sorted((str(2*int(locus[0])-1),str(2*int(locus[0])))))
            loc_orient = 0
            if (locus in block_indices and tuple(locus[::-1]) in block_indices):
                loc_orient = 0
            elif locus in block_indices:
                loc_orient = 1
            elif tuple(locus[::-1]) in block_indices:
                loc_orient = -1
            else:
                loc_orient = 0
            in_locus = []
            in_species = {}
            list_in_species = {}
            adj_in_species = []
            possible_content = []
            local_copies = {}
            kept_species = []
            addable = False
        elif not(line.isspace()):
            sp,chr = line[:line.find(':')].split('.')
            if sp not in in_species:
                in_species[sp] = {}
            if chr not in in_species[sp]:
                in_species[sp][chr] = set([])
            if sp not in list_in_species:
                list_in_species[sp] = {}
            if chr not in list_in_species[sp]:
                list_in_species[sp][chr] = []
            temp_content = line[line.find(':')+1:].strip().split(' ')[:-1]
            list_in_species[sp][chr] += temp_content
            orient = line[line.find(':')+1:].strip().split(' ')[-1]
            if locus not in big_adj_orients:
                big_adj_orients[locus] = {}
            if sp not in big_adj_orients[locus]:
                big_adj_orients[locus][sp] = {}
            if chr not in big_adj_orients[locus][sp]:
                if orient == '+':
                    big_adj_orients[locus][sp][chr] = 1
                elif orient == '-':
                    big_adj_orients[locus][sp][chr] = -1
            if temp_content == ['']:
                temp_content.pop()
            kept_species.append(sp)
            addable = True                
            for k in range(len(temp_content)):
                if temp_content[k] not in local_copies:
                    local_copies[temp_content[k]] = {}
                    local_copies[temp_content[k]][sp] = 0
                if sp not in local_copies[temp_content[k]]:
                    local_copies[temp_content[k]][sp] = 0
                local_copies[temp_content[k]][sp] += 1
                if k%2 == 1 and k < len(temp_content)-1:
                    in_species[sp][chr].add(tuple(sorted((temp_content[k],temp_content[k+1]))))
                    if tuple(sorted((temp_content[k],temp_content[k+1]))) not in possible_content:
                        possible_content.append(tuple(sorted((temp_content[k],temp_content[k+1]))))
                    if temp_content[k] < temp_content[k+1]:
                        if (temp_content[k],temp_content[k+1]) not in adj_orients:
                            adj_orients[(temp_content[k],temp_content[k+1])] = {}
                            adj_orients[(temp_content[k],temp_content[k+1])][-1] = 0
                            adj_orients[(temp_content[k],temp_content[k+1])][1] = 0
                        if orient == '-':
                            adj_orients[(temp_content[k],temp_content[k+1])][1] += 1
                        elif orient == '+':
                            adj_orients[(temp_content[k],temp_content[k+1])][-1] += 1
                    if temp_content[k] > temp_content[k+1]:
                        if (temp_content[k],temp_content[k+1]) not in adj_orients:
                            adj_orients[(temp_content[k],temp_content[k+1])] = {}
                            adj_orients[(temp_content[k],temp_content[k+1])][-1] = 0
                            adj_orients[(temp_content[k],temp_content[k+1])][1] = 0
                        if orient == '-':
                            adj_orients[(temp_content[k],temp_content[k+1])][-1] += 1
                        elif orient == '+':
                            adj_orients[(temp_content[k],temp_content[k+1])][1] += 1
                    elif temp_content[k] == temp_content[k+1]:
                        if (temp_content[k],temp_content[k+1]) not in adj_orients:
                            adj_orients[(temp_content[k],temp_content[k+1])] = {}
                            adj_orients[(temp_content[k],temp_content[k+1])][-1] = 0
                            adj_orients[(temp_content[k],temp_content[k+1])][1] = 0

if len(in_species) > 0 and addable:
    if locus not in content:
        content[locus] = findParsimoniousContent(locus,tree_file,in_species,possible_content,local_copies,copynumbers)
        unparsed_content[locus] = list_in_species 
    else:
        temp_content = findParsimoniousContent(locus,tree_file,in_species,possible_content,local_copies,copynumbers)
        unparsed_content[locus] = list_in_species 
        content[locus] += temp_content[0]
        for k in temp_content[1]:
            if k in content[locus][1] and \
                    content[locus][1][k] < temp_content[1][k]:
                content[locus][1][k] = temp_content[1][k]
            elif k not in content[locus][1]:
                content[locus][1][k] = temp_content[1][k]
        for m in temp_content[2]:
            if m not in content[locus][2]:
                content[locus][2] = temp_content[2][m]
            else:
                for s in temp_content[2][m]:
                    if s not in content[locus][2][m]:
                        content[locus][2][m].append(m)


marker_assoc = {}
marker_list = {}
marker_map = {}
current_localization = {}
adjGraph = nx.Graph()
gene_to_loci = {}
for c in content:
    if c in block_indices:
        internal_labelling = augmentGraph(adjGraph,marker_list,\
                content[c][0],content[c][1],\
                copynumbers,current_localization)
        marker_assoc[c] = internal_labelling
        car_id = block_indices[c].keys()[0]
        for g in internal_labelling:
            gl = internal_labelling[g]
            addToGeneLocus(gene_to_loci,gl,car_id,c)
            addToGeneLocus(gene_to_loci,makeMate(gl),car_id,c)

adjs = {}
with open(sys.argv[3]) as adj_file:
    for line in adj_file:
        weight = float(line[:line.find(';')].split('|')[-1].strip())
        temp_content = line[line.find(':')+1:line.find('#')].strip().split(' ')
        adjs[(temp_content[0],temp_content[1])] = weight#*.1
        adjs[(temp_content[1],temp_content[0])] = weight#*.1


segment_size = int(sys.argv[9])
window_size = int(sys.argv[10])
covered = set([])
remove_set = []
empty_set = {}
added_adjs = set([])
merged_to = {}

for current_car in car_dict:            
    i = 1
    car = car_dict[current_car]
    addExtraContent(adjGraph,current_car,car,unparsed_content,content,tree_file,marker_assoc,gene_to_loci,\
            window_size,segment_size,copynumbers,big_adj_orients,current_localization,adjs)
    while i < len(car)-1:
        if tuple(sorted((car[i-1],car[i]))) in content:
            c0 = tuple(sorted((car[i-1],car[i])))
        else:
            i += 1
            continue
        prev = i
        while i < len(car)-1 and \
                tuple(sorted((car[i],car[i+1]))) not in content:
                i += 1
        j = 0 
        window = []
        while i+j+1 < len(car) and j < window_size:
            if i+j+1 < len(car) and \
                    tuple(sorted((car[i+j+1],car[i+j]))) in content:
                window.append(tuple(sorted((car[i+j],car[i+j+1]))))
            j += 1            
        for c1 in window:    
            remove_set_temp, used_temp = addAdjs(adjGraph,c0,c1,adjs,marker_assoc,gene_to_loci,\
                    copynumbers,covered,merged_to)
            remove_set += remove_set_temp
            added_adjs = added_adjs.union(used_temp)
        sys.stdout.flush()
        if i == prev:
            i += 1
        else:
            i = prev+1

for r in remove_set:
    if r in adjGraph:
        adjGraph.remove_node(r)
        del gene_to_loci[r]
    if makeMate(r) in adjGraph:    
        adjGraph.remove_node(makeMate(r))
        del gene_to_loci[makeMate(r)]

profiles = {}
with open(sys.argv[5]) as pro_file:
    for line in pro_file:
        if '>' in line:
            id,num = map(int,line[1:].strip().split(' '))
            t = str(2*id-1)
            h = str(2*id)
            profiles[t] = num
            profiles[h] = num

accounted_for = set([])
for v in adjGraph:
    nbs = adjGraph[v].keys()
    nbs.sort(key=lambda x: x[0])
    deleteRedundantEdges(adjGraph,v,nbs,gene_to_loci,block_indices)



for e in adjGraph.edges():
    accounted_for.add(tuple(sorted([e[0][0],e[1][0]])))

for a in adjs:
    if tuple(sorted(a)) not in accounted_for:#covered and a not in added_adjs:
        v0,v1 = a
        if tuple(sorted([a])) not in adj_orients:
            adj_orients[tuple(sorted([a]))] = {}
            adj_orients[tuple(sorted([a]))][1] = 0
            adj_orients[tuple(sorted([a]))][-1] = 0
        i = 0
        if v0 not in [x[0] for x in adjGraph]:
            current_localization[v0] = 0
            current_localization[makeMate(v0)] = 0
            adjGraph.add_node((v0,current_localization[v0]))
            adjGraph.add_node(makeMate((v0,current_localization[v0])))
            copynumbers[makeMate((v0,0))] = profiles[v0]
            copynumbers[(v0,0)] = profiles[v0]
        if v1 not in [x[0] for x in adjGraph]:
            current_localization[v1] = 0
            current_localization[makeMate(v1)] = 0
            adjGraph.add_node((v1,current_localization[v1]))
            adjGraph.add_node(makeMate((v1,current_localization[v1])))
            copynumbers[(v1,0)] = profiles[v1]
            copynumbers[makeMate((v1,0))] = profiles[v1]
        #################################################################            
        while i <= current_localization[v0]:
            if (v0,i) in adjGraph.nodes():
                if copynumbers[(v0,i)] > profiles[v0]:
                    copynumbers[makeMate((v0,i))] = profiles[v0]
                j = 0
                while j <= current_localization[v1]:
                    if (v1,j) in adjGraph.nodes():
                        if copynumbers[(v1,j)] > profiles[v1]:
                            copynumbers[(v1,j)] = profiles[v1]
                            copynumbers[makeMate((v1,j))] = profiles[v1]
                        adjGraph.add_edge((v0,i),(v1,j),weight=adjs[a])
                    j += 1
            i += 1
        #################################################################            

orig_nodes = adjGraph.nodes()
for v in orig_nodes:
    if v not in gene_to_loci:
        addToWindow(adjGraph,v,marker_assoc,gene_to_loci,block_indices)

for v in profiles:
    i = 0
    v_content = []
    if v not in current_localization:
#        print "Err",v
        continue
    while i <= current_localization[v]:
        if (v,i) in adjGraph.nodes():
            v_content.append(i)
        i += 1            
    if profiles[v] < sum([copynumbers[(v,j)] for j in v_content]):
        balanceCopies(adjGraph,v,v_content,copynumbers,profiles[v],\
                marker_assoc,gene_to_loci)


o = open(sys.argv[6],'w')
done = []
to_add = []
for c in copynumbers:
    if makeMate(c) not in copynumbers:
        to_add.append(makeMate(c))

for c in to_add:
    copynumbers[c] = copynumbers[makeMate(c)]    
del to_add

for c in copynumbers:
    if c not in adjGraph:
        continue
    undoubled = int(c[0])
    if undoubled%2 != 0:
        undoubled += 1
    undoubled = undoubled >> 1
    if (undoubled,c[1]) not in done:
        o.write('>'+str(undoubled)+'.'+str(c[1])+' '+str(copynumbers[c])+'\n')
        done.append((undoubled,c[1]))

o.close()
car_tiling = findSubCARs(adjGraph,marker_assoc,gene_to_loci,\
        car_dict, copynumbers, int(sys.argv[9]),adj_orients)


o = open(sys.argv[7],'w')
empties = 0
for c in car_dict:
    o.write('#'+str(c)+'\n')
    i = 0
    for paths in car_tiling[c]:
        o.write('Tiling: '+' '.join(car_dict[c][i:i+int(sys.argv[9])])+'\n')
        if len(paths) > 0:
            o.write('\n'.join('\t'+' '.join(['.'.join(map(str,x)) for x in p])for p in paths if len(p) > 0 )+'\n')
        else:            
            empties += 1
        o.write('$$\n')
        i += 1
    o.write('\n')
#o.write(str(empties))
o.close()
o = open(sys.argv[8],'w')

for c in marker_assoc:
    o.write('>'+' '.join(c)+'\n')
    o.write('\n'.join([marker_assoc[c][u][0]+'.'+str(marker_assoc[c][u][1]) for u in marker_assoc[c]])+'\n')
o.close()    

