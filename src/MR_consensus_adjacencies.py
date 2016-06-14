#!/usr/bin/env

import sys
from marker_string import *
from adjtools import *
import networkx as nx
import linearization as lz
import random
import kmp
import copy
from  graph_traverser import *

def getAllPairwiseAdjs(marker_list):
    unfiltered_putative_adjs = []
    for i in range(len(marker_list)):
        m = marker_list[i]
        m_id,m_loc = m.split('.')
        m_loc = int(m_loc)
        m_tail,m_head = readUndoubledMarker(m_id,m_loc)
        for j in range(i,len(marker_list)):
            mm = marker_list[j]
            mm_id,mm_loc = mm.split('.')
            mm_loc = int(mm_loc)
            mm_tail,mm_head = readUndoubledMarker(mm_id,mm_loc)
            for l in [m_tail,m_head]:
                for k in [mm_tail,mm_head]:
                    if tuple(sorted([l,k])) not in unfiltered_putative_adjs:
                        unfiltered_putative_adjs.append(tuple(sorted((l,k))))
    return unfiltered_putative_adjs                        

def getBlockAdjs(m0,m1,assoc_dict):
    current_block = tuple(sorted((m0,m1))) 
    if current_block in assoc_dict:
        return getAllPairwiseAdjs(assoc_dict[current_block])
    else:
        return KeyError

# path: list of gene extremities in order
# local_mults: empty dict to be populated
# path_adjs: list of inter-gene adjacencies in path, each in lex order
def parsePath(path,local_mults):
    path_adjs = []
    for i in range(1,len(path)-2,2):
        if path[i] not in local_mults:
            local_mults[path[i]] = 0
        if makeMate(path[i]) not in local_mults:
            local_mults[makeMate(path[i])] = 0
        local_mults[path[i]] += 1
        local_mults[makeMate(path[i])] += 1
        path_adjs.append(tuple(sorted((path[i],path[i+1]))))
    if path[-1] not in local_mults:
        local_mults[path[-1]] = 0
    if makeMate(path[-1]) not in local_mults:
        local_mults[makeMate(path[-1])] = 0
    local_mults[path[-1]] += 1
    local_mults[makeMate(path[-1])] += 1
    return path_adjs

# window_paths: list of paths in a window
# start_index: index of window
# all_mults: dict to be populated with estimated copy numbers, 
#            given as max of local copy numbers; updated by window
# window_adjs: list of adjacencies in the window
def parseWindow(window_paths,start_index,all_mults,marker_assoc,gene_assoc):
    local_mults = {}
    window_adjs = []
    for p in window_paths:
        window_adjs += parsePath(p,local_mults)
    for k in local_mults:
        if k in all_mults and all_mults[k] < local_mults[k]:
            all_mults[k] = local_mults[k]
        elif k not in all_mults:
            all_mults[k] = local_mults[k]
        if makeMate(k) in all_mults and all_mults[makeMate(k)] < local_mults[makeMate(k)]:
            all_mults[makeMate(k)] = local_mults[makeMate(k)]
        elif makeMate(k) not in all_mults:
            all_mults[makeMate(k)] = local_mults[makeMate(k)]
    return window_adjs

def allDoubled(marker_list):
    to_append = []
    for m in marker_list:
        if '.' in m:
            m_id, loc = m.split('.')
            loc = int(loc)
            to_append.append(readDoubledMarker(m))
        else:
            m_id = m
            to_append.append(readDoubledMarker(m))
    return to_append

def mustDelete(G,u,nbhood,c,already_in,already_counted):
    nbhood.sort(key = lambda x: G[u][x]['weight'])
    deletable = []
    constraint = c
    if makeMate(u) in nbhood:
        constraint = c+1
    for v in nbhood:
        if G[u][v]['color'] == 'r': 
            if (u,v) not in already_in and (v,u) not in already_in and \
                    len(nbhood)-len(deletable) > constraint:
                deletable.append(v)
                if v in already_counted:
                    already_counted[v] += 1
                else:
                    already_counted[v] = 1
            elif len(nbhood)-len(deletable) <= constraint:
                break
    if u not in already_counted:
        already_counted[u] = 0
    already_counted[u] += len(deletable)
    return [(u,v) for v in deletable]



def deletableEdges(G,u,mults,already_in,already_counted):
    c = mults[u]
    nbhood = G[u].keys()
    already_deleted = 0
    if u in already_counted:
        already_deleted = already_counted[u]
    if len(nbhood)-already_deleted <= c:
        return []
    else:
        return mustDelete(G,u,nbhood,c+already_deleted,already_in,already_counted)

def isDeletable(G,e,mults,deletion_list):
    u,v,k = e
    u_m = makeMate(u)
    v_m = makeMate(v)
    if (redLen(G,u)-deletion_list[u] > mults[u] or \
            redLen(G,v)-deletion_list[v] > mults[v]):
        return True
    return False

def reduceMults(e,deletion_list):
    for v in e[:2]:
        deletion_list[v] += 1

def deletableEdgesParent(G,mults):
    temp_mults = mults
    for x in G:
        if x not in temp_mults:
            temp_mults[x] = 1
            temp_mults[makeMate(x)] = 1
    to_delete = []
    edge_list = sorted([x for x in G.edges(keys=True) if G[x[0]][x[1]][x[2]]['color']=='r'],\
            key=lambda x:G[x[0]][x[1]][x[2]]['weight'],reverse=True)
#    print "Checking deletion",[G[e[0]][e[1]][e[2]]['weight'] for e in edge_list]
    deletion_list = dict([(x,0) for x in G])
    #print mults[('641',14)]
    #print deletion_list[('641',14)]
    while any([redLen(G,x)-deletion_list[x] > mults[x] for x in G]):
        e = edge_list.pop()
        if isDeletable(G,e,mults,deletion_list):
            if G[e[0]][e[1]][e[2]]['color'] != 'b':
#                print "Deleting",e,G[e[0]][e[1]][e[2]]['weight']
                to_delete.append(e)
                reduceMults(e,deletion_list)
    return to_delete

def findCommonWindows(u,v,gene_associations,car_id):
#    print u,v, gene_associations[u],gene_associations[v],car_id
    if car_id not in gene_associations[u] or car_id not in gene_associations[v]:
        return sys.maxint
    len_common_wins = len(\
            set(\
            range(min(gene_associations[u][car_id]),\
            max(gene_associations[u][car_id])+1)\
            )\
            &\
            set(\
            range(min(gene_associations[v][car_id]),\
            max(gene_associations[v][car_id])+1)\
            )\
            )
    if len_common_wins == 0:
        len_common_wins = sys.maxint
#                print u,v,gene_associations[u][car_id],gene_associations[v][car_id] 
    return len_common_wins



def addMateEdge(G,v):
    to_add = True
    for u in G[v]:
        if u != makeMate(v):
            continue
        for k in G[v][u]:
            if G[v][u][k]['color'] == 'b':
                to_add = False
                break
        if not(to_add):
            break
    if to_add:        
        G.add_edge(v,makeMate(v),weight=sys.maxint,color='b')


def repath(path_seq,part_seq):
    i = 0
    j = 0
    new_path_seq = []
    prev_reverse = False
    curr_reverse = False
    added = False
    prev_pos = -1
    printed_0 = False
    while i < len(part_seq)-1 and j < len(path_seq):
        part0 = path_seq[j:part_seq[i]]
        part1 = path_seq[part_seq[i]+1:part_seq[i+1]]
        informative = False
        if not(added) and part0[0] in part1:
            prev_reverse = True
        if not(added) and prev_reverse:
            new_path_seq += part0[::-1] + ['|']
        elif not(added):
            new_path_seq += part0 + ['|']
        if part1[-1] in part0:
            informative = True
            curr_reverse = True
        elif part1[0] in part0:
            curr_reverse = False
            informative = True
        if informative and curr_reverse:
            new_path_seq += part1[::-1] +['|']
            added = True
        elif informative and not(curr_reverse):
            new_path_seq += part1 +['|']
            added = True
        else:
            added = False
        j = part_seq[i]+1 
        i += 1
    if i > 0:
        part0 = path_seq[part_seq[i-1]+1:part_seq[i]]
        part1 = path_seq[part_seq[i]+1:]
    else:
        part0 = path_seq[:part_seq[i]]
        part1 = path_seq[part_seq[i]+1:]
    informative = False
    if not(added) and len(part0) > 0 and part0[0] in part1:
        prev_reverse = True
    if not(added) and prev_reverse:
        new_path_seq += part0[::-1] + ['|']
    elif not(added):
        new_path_seq += part0 + ['|']
    if part1 != []:    
        if part1[-1] in part0:
            informative = True
            curr_reverse = True
        elif part1[0] in part0:
            curr_reverse = False
            informative = True
        if informative and curr_reverse:
            new_path_seq += part1[::-1]
#        elif informative and not(curr_reverse):
#            new_path_seq += part1
        else:
            new_path_seq += part1 
    return new_path_seq 

def weightComponent(contents,g):
    try:
        wt = len(set(contents)&g)
        return wt
    except:
        return 0

# Given a block (m0,m1) in the coarse CAR,
# check how many markers in the current component are common to the 
# predicted block content.
# Higher weight 
def weightComponents(m0,m1,current_index,marker_dict,G_ccs,comp_order):
    try:
        contents = marker_dict[tuple(sorted([m0,m1]))]
    except KeyError:
        contents = []
    if len(contents) > 0:
        for g in G_ccs:
            wt = weightComponent(contents,set(g.nodes()))
            if wt > 0 and (wt >= len(contents)>>2 or wt == len(g.nodes())):
                if g not in comp_order:
                    comp_order[g] = []
                comp_order[g].append(current_index)    

def addTandems(adj_list):
    i = 0
    while i < len(adj_list)-1:
        e0,e1 = adj_list[i]
        if (e0 in adj_list[i+1] and makeMate(e0) not in adj_list[i+1]) and \
            makeMate(e1) not in adj_list[i+1]:
            adj_list.insert(i+1,(makeMate(e0),e0))
        elif (e1 in adj_list[i+1] and makeMate(e1) not in adj_list[i+1]) and \
                makeMate(e0) not in adj_list[i+1]:
            adj_list.insert(i+1,(makeMate(e1),e1))
        i += 1            



#for i in range(2,2+n_res):
#    car_paths,window_list,win_size = extractConstructedPaths(sys.argv[i])
#    res_paths[win_size] = car_paths
#
def extractConstructedPaths(file_name):
    car_paths = {}
    car_window_list = {}
    car_id = ''
    path_list = []
    win_size = -1
    with open(file_name) as constructed_paths_file:
        for line in constructed_paths_file:
            if not(line.isspace()):
                if line[0] == '#':
                    car_id = line[1:].strip()
                    car_paths[car_id] = []
                    car_window_list[car_id] = []
                elif line[0] == 'T':
                    car_window_list[car_id].append(line[line.find(':')+1:].strip().split(' '))
                    if len(line[line.find(':')+1:].strip().split(' ')) > win_size:
                        win_size =len( line[line.find(':')+1:].strip().split(' '))
                    window_paths = []
                elif line[0] == '$':
                    car_paths[car_id].append(window_paths)
                elif line[0] == '\t':
                    window_paths.append(map(readDoubledMarker,line.strip().split(' ')))
                else:
                    continue
    return car_paths,car_window_list,win_size





def traverseGraph(G,adj_windows,window_adjs):
    start_vertex = None
    k = 0
    for e in window_adjs[k]:
        dothings()


def determineOrder(adjs,windows):
    for i in range(len(adjs)-1):
        a = adjs[i]
        j = i+1
        while j < len(adjs) and \
                (makeMate(a[0]) not in adjs[j] and \
                makeMate(a[1]) not in adjs[j]):
            j += 1
        if j < len(adjs):
            adjs[i+1:j+1] = adjs[i+1:j+1][::-1]
            adjs[i+2:j+1] = adjs[i+2:j+1][::-1]
#            adjs[i+1],adjs[j] = adjs[j],adjs[i+1]

def resort(current_slice, window_dict):
    for i in range(1,len(current_slice)):
        e = current_slice[i-1]
        j = i
        while j < len(current_slice) and \
            (makeMate(e[0]) not in current_slice[j] or \
            makeMate(e[1]) not in current_slice[j]):
                j += 1
        if j < len(current_slice):
            current_slice[i],current_slice[j] = current_slice[j],current_slice[i]
    return current_slice

def populateStack(nb_list,stack,add_list,gene_assoc,car_id,to_check):
    for x in nb_list:
        if x not in gene_assoc and x not in stack:
            stack.append(x)
            stack.append(makeMate(x))
            to_check.append(x)
            to_check.append(makeMate(x))
        elif x in gene_assoc and \
                car_id not in gene_assoc[x] and \
                x not in stack:
            stack.append(x)
            stack.append(makeMate(x))
            to_check.append(x)
            to_check.append(makeMate(x))
        elif x in gene_assoc and car_id in gene_assoc[x]:
            add_list += gene_assoc[x][car_id]
    
def populateAssocDictInner(stack,gene_assoc,to_add,car_id):
    for x in stack:
        if x not in gene_assoc:
            gene_assoc[x] = {}
            gene_assoc[makeMate(x)] = {}
        if car_id not in gene_assoc:
            gene_assoc[x][car_id] = []
            gene_assoc[makeMate(x)][car_id] = []
        gene_assoc[x][car_id] += to_add
        gene_assoc[makeMate(x)][car_id] += to_add
        gene_assoc[x][car_id] = list(set(gene_assoc[x][car_id]))
        gene_assoc[makeMate(x)][car_id] = list(set(gene_assoc[makeMate(x)][car_id]))

def populateAssocDict(G,v,gene_assoc,car_id):
    add_to_left = []
    add_to_right = []
    left_stack = [v]
    right_stack = [makeMate(v)]
    left_check = [v]
    right_check = [makeMate(v)]
    while left_check != [] or right_check != []:
        if left_check != []:
            l = left_check.pop()
            nbs_left = [u for u in G[l] for k in G[l][u] if (G[l][u][k]['color'] == 'r' and u != makeMate(l))] 
            populateStack(nbs_left,left_stack,add_to_left,gene_assoc,car_id,left_check)
        if right_check != []:
            r = right_check.pop()
            nbs_right = [u for u in G[r] for k in G[r][u] if (G[r][u][k]['color'] == 'r' and u != v)]
            populateStack(nbs_right,right_stack,add_to_right,gene_assoc,car_id,right_check)
    populateAssocDictInner(left_stack,gene_assoc,add_to_left,car_id)    
    populateAssocDictInner(right_stack,gene_assoc,add_to_right,car_id)    

def isAllMatesIn(test):
    all_in = True
    for t in test:
        t0,t1 = t
        all_in = all_in and (makeMate(t0))
        all_in = all_in and (makeMate(t1))
        if not(all_in):
            return False
    return True        

def insertionSort(adjacency_list,window_dict):
    new_list = []
    i = 0
    while i < len(adjacency_list):
        j = i
        reference_list = window_dict[adjacency_list[i]]
        st_win = reference_list[0]
        end_win = reference_list[-1]
        while j < len(adjacency_list) and \
                ((window_dict[adjacency_list[j]][0] >= st_win and \
                window_dict[adjacency_list[j]][-1] <= end_win) or \
                not(isAllMatesIn(adjacency_list[i:j]))):
            j += 1                    
        slice = adjacency_list[i:j]
        sorted_slice = resort(adjacency_list[i:j],window_dict)
        new_list += sorted_slice
        i = j
    return new_list
        
# car_paths: dict with car ids as keys, list of list of paths as values (in window order, then arbitrary)
# car_window_list: dict with car ids as keys, list of windows as values
# win_size: positive integer

car_paths, car_window_list, win_size = extractConstructedPaths(sys.argv[3])

# marker_associations: dict with locus id (lex order) as keys, list of gene extremities as values

marker_associations = {}
m_id = ''
contents = []
with open(sys.argv[1]) as association_file:
    for line in association_file:
        if line.isspace():
            continue
        if '>' in line:
            if m_id != '':
                marker_associations[m_id] = allDoubled(contents)
            m_id = tuple(sorted(line[1:].strip().split()))
            contents = []
        else:
            contents.append(line.strip())
marker_associations[m_id] = allDoubled(contents)

# cars: dict with car ids as keys, list of doubled synteny loci as values
cars = {}
gene_associations = {}
all_loci = set([])
with open(sys.argv[2]) as car_file:
    car_id = ''
    for line in car_file:
        if not(line.isspace()):
            if line[0] == '#':
                car_id = line[1:].strip()
            elif line[0] == '_':
                cars[car_id] = line.strip().split(' ')[1:-1]
#                print line.strip().split(' ')
#                print cars[car_id]
                current_car = cars[car_id]
                for i in range(len(current_car)-1):
                    locus = tuple(sorted([current_car[i],current_car[i+1]]))
                    all_loci.add(locus)
                    if locus in marker_associations:
                        for g in marker_associations[locus]:
                            if g not in gene_associations:
                                gene_associations[g] = {}
                                gene_associations[makeMate(g)] = {}
                            if car_id not in gene_associations[g]:
                                gene_associations[g][car_id] = []
                                gene_associations[makeMate(g)][car_id] = []
                            if i not in gene_associations[g][car_id]:
                                gene_associations[g][car_id].append(i)
                            if i not in gene_associations[makeMate(g)][car_id]:
                                gene_associations[makeMate(g)][car_id].append(i)

#window_length = int(sys.argv[6])

to_remove = []
for l in marker_associations:
    if l not in all_loci:
        to_remove.append(l)
for l in to_remove:
    del marker_associations[l]

for g in gene_associations:
    for car_id in gene_associations[g]:
        interval = sorted(gene_associations[g][car_id])
        gene_associations[g][car_id] = interval#[min(interval), max(interval)]

copy_numbers = {}
with open(sys.argv[4]) as cp_file:
    for line in cp_file:
        if not(line.isspace()):
            m,num = line[1:].strip().split(' ')
            m,loc = m.split('.')
            loc = int(loc)
            for mm in readUndoubledMarker(m,loc):
                copy_numbers[mm] = int(num)

segment_length = int(sys.argv[5])
window_length = int(sys.argv[6])

for g in gene_associations:
    for car in gene_associations[g]:
        if len(gene_associations[g][car]) < 3:
            continue
        end0 = min(gene_associations[g][car])
        end1 = max(gene_associations[g][car])
        while abs(end0-end1) > window_length:
            if abs(end0-gene_associations[g][car][1]) > abs(end1-gene_associations[g][car][-2]):
                gene_associations[g][car] = gene_associations[g][car][1:]
                gene_associations[makeMate(g)][car] = gene_associations[makeMate(g)][car][1:]
                locus = tuple(sorted([cars[car][end0],cars[car][end0+1]]))
                if g in marker_associations[locus]:
                    marker_associations[locus].remove(g)
                if makeMate(g) in marker_associations[locus]:
                    marker_associations[locus].remove(makeMate(g))
            else:# abs(end0-gene_associations[g][car][1]) > abs(end1-gene_associations[g][car][-2]):
                gene_associations[g][car] = gene_associations[g][car][:-1]
                gene_associations[makeMate(g)][car] = gene_associations[makeMate(g)][car][:-1]
                locus = tuple(sorted([cars[car][end1],cars[car][end1+1]]))
                if g in marker_associations[locus]:
                    marker_associations[locus].remove(g)
                if makeMate(g) in marker_associations[locus]:
                    marker_associations[locus].remove(makeMate(g))
            end0 = min(gene_associations[g][car])
            end1 = max(gene_associations[g][car])

# adj_windows: dict with lex ordered adjs as keys 
for c in car_window_list:
    adj_windows = {}
    window_adjs = [] 
    new_mults = {}
    marker_windows = {} # Associates a marker to a list of windows
    for i in range(0,len(car_window_list[c])):
        paths_list = car_paths[c][i]
        window_adjs.append(parseWindow(paths_list,i,\
                new_mults,marker_associations,gene_associations))
        for a in window_adjs[-1]:
            if a not in adj_windows:
                adj_windows[a] = []
            adj_windows[a].append(i)
            if a[0] not in marker_windows:
                marker_windows[makeMate(a[0])] = []
                marker_windows[a[0]] = []
            marker_windows[a[0]].append(i)
            marker_windows[makeMate(a[0])].append(i)
            if a[1] not in marker_windows:
                marker_windows[makeMate(a[1])] = []
                marker_windows[a[1]] = []
            marker_windows[makeMate(a[1])].append(i)
            marker_windows[a[1]].append(i)
    geneGraph = nx.MultiGraph()
    for a in adj_windows:
        common_windows = findCommonWindows(a[0],a[1],gene_associations,c)
        geneGraph.add_edge(a[0],a[1], weight=len(adj_windows[a])*1./common_windows, color='r')
        addMateEdge(geneGraph,a[0])
        addMateEdge(geneGraph,a[1])
    nodes_to_remove = []        
    for v in geneGraph:
        if v not in gene_associations:
            gene_associations[v] = {}
            gene_associations[makeMate(v)] = {}
        if c not in gene_associations[v]:
            nodes_to_remove.append(v)
            nodes_to_remove.append(makeMate(v))
    geneGraph.remove_nodes_from(nodes_to_remove)            
    # Using calculated localized mults        
    to_delete = deletableEdgesParent(geneGraph,copy_numbers)
    sys.stdout.flush()
    geneGraph.remove_edges_from(to_delete)
    geneGraph_cc = []
    #print nx.number_connected_components(geneGraph)
    for g in nx.connected_component_subgraphs(geneGraph):
        geneGraph_cc.append(g)
    for q in [1]:
        print '#'+c
        if len(geneGraph.nodes()) == 0:
            continue
        path = componentTraverser(geneGraph_cc,c,gene_associations,copy_numbers)
        print ' '.join(['.'.join(map(str,x)) for x in path])+" $"
    
