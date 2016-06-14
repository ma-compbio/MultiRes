#!/usr/bin/env
import sys
from marker_string import *
import copy 
from DataStructures.binary_heap import *
import networkx as nx

def findEndOrder(end_list,gene_assoc,car_id):
    end_win = {}
    max_len = 0
    for v in end_list:
        if car_id in gene_assoc[v]:
            end_win[v] = sorted(gene_assoc[v][car_id])
            gene_assoc[v][car_id].sort()
            if len(end_win[v]) > max_len:
                max_len = len(end_win[v]) 
    sorted_ends = sorted(end_win.keys(), \
            key=lambda x: min(end_win[x]))
    return end_win,sorted_ends


def findEnds(G,gene_assoc,mults,car_id):
    ends = []
    if len(G) == 2:
        ends = G.nodes()
        return ends
    for v in G:
        if G.degree(v) == 1 and mults[v] < 2:
            ends.append(v)            
        elif makeMate(v) in [u for u in G[v] \
                for k in G[v][u] if G[v][u][k]['color'] == 'r']:
            red_non_tandem = [u for u in G[v] \
                for k in G[v][u] if \
                (G[v][u][k]['color'] == 'r' and u != makeMate(v))]
            if len(red_non_tandem) == 1 and mults[red_non_tandem[0]] < 2:
                ends.append(v)                
    if ends == []:
        to_delete = None
        max_discord = 0
        mwe = None
        min_wt = sys.maxint
        for v in G:
            for u in G[v]:
                for k in G[v][u]:
                    if G[v][u][k]['color'] == 'r':
                        if abs(gene_assoc[v][car_id][0]-gene_assoc[u][car_id][0]) \
                                > max_discord:
                            max_discord = abs(gene_assoc[v][car_id][0]-gene_assoc[u][car_id][0])
                            to_delete = (v,u,k)
                        if G[v][u][k]['weight'] < min_wt:
                            min_wt = G[v][u][k]['weight']
                            mwe = (v,u,k)
        if to_delete and max_discord > 5:
            G.remove_edge(to_delete[0],to_delete[1],to_delete[2])
            ends.append(to_delete[0])
            ends.append(to_delete[1])
        else:
            ends.append(mwe[0])
            ends.append(mwe[1])
            G.remove_edge(mwe[0],mwe[1],mwe[2])
        all_ends = []
    end_windows,sorted_ends = findEndOrder(ends,gene_assoc,car_id)
    true_ends = []
    for e in sorted_ends:
        if makeMate(e) not in end_windows and \
                G.degree(makeMate(e)) <= 2:
            true_ends.append(e)
    if true_ends == []:
        true_ends = ends
    true_ends.sort(key=lambda x: end_windows[x][0])
    return true_ends

def rsiLocater(G,path_list,mults,rsi_dict={}):
    for p in path_list:
        i = 0
        while i < len(p)-1:
            start_rsi = False
            if mults[p[i]] == 1 and \
                    mults[p[i+1]] > 1:
                start_rsi = True                        
            if start_rsi:
                j = i+1
                while j < len(p) or \
                        (j < len(p)-1 and mults[p[j+1]] > 1):
                    j += 1                            
                if j < len(p):
                    rsi = tuple(sorted(p[i:j+1]))
                    if rsi not in rsi_dict:
                        rsi_dict[rsi] = 0
                    rsi_dict[rsi] += 1


def extendCluster(G,v,mults):
    to_traverse = [v]
    traversed = set([v])
    while to_traverse != []:
        current_vertex = to_traverse.pop()
        for u in G[current_vertex]:
#            print "NBS",current_vertex,u
            if mults[u] > 1 and \
                    u not in traversed:
                to_traverse.append(u)
                traversed.add(u)
            if mults[u] > 1 and \
                    makeMate(u) not in traversed:
                to_traverse.append(makeMate(u))
                traversed.add(makeMate(u))
    return traversed 

def repeatClusterLocation(G,marker_windows,mults):
    in_cluster = set([])
    clusters = []
    for v in G:
        if mults[v] > 1 and \
                v not in in_cluster and \
                makeMate(v) not in in_cluster:
            R = extendCluster(G,v,mults)
            #print R
            in_cluster = in_cluster | R
            #print in_cluster
            clusters.append(R)
    return clusters            

def isCloser(expect,gene_assoc,car_id,m0_list,m1_list):
    min_approach0 = minApproach(m0_list,gene_assoc,car_id,expect)
    min_approach1 = minApproach(m1_list,gene_assoc,car_id,expect)
    if min_approach0 < min_approach1:
        return True
    elif min_approach0 < min_approach1:
        return False
    return 0 

def minApproach(nb_list,gene_assoc,car_id,expect):
    minimum_app = sys.maxint
    for n in nb_list:
        new_l = gene_assoc[n][car_id]
#        print "Searching",n,new_l
        for x in new_l:
            for e in expect:
                if x > e and x-e < minimum_app:
                    minimum_app = x-e
    return minimum_app

#    end_list = traverser(G,marker_windows,mults)
def nb_search(G,gene_assoc,m0,m1,car_id,prospect,signal):
    expect = gene_assoc[prospect][car_id]
    m0_list = redNbhood(G,m0) + redNbhood(G,makeMate(m0))
    m1_list = redNbhood(G,m1) + redNbhood(G,makeMate(m1))
    is_closer = 0 
    used = []
    while is_closer == 0 and m0_list != [] and m1_list != []:
#        print "Test",m0,m1,m0_list,m1_list
        if prospect in m0_list:
            m0_list.remove(prospect)
        if prospect in m1_list:
            m1_list.remove(prospect)
        is_closer = isCloser(expect,gene_assoc,car_id,m0_list,m1_list)
        if is_closer != 0:
            break
        else:
#            print "Search tree needed"
            tm0_list = m0_list[:]
            tm1_list = m1_list[:]
#            print tm0_list,m0_list
#            print tm1_list,m1_list
            m0_list = []
            m1_list = []
            for x in tm0_list:
                if makeMate(x) not in used:
                    #print redNbhood(G,makeMate(x))
                    m0_list += redNbhood(G,makeMate(x))
            for x in tm1_list:
                if makeMate(x) not in used:
                    #print redNbhood(G,makeMate(x))
                    m1_list += redNbhood(G,makeMate(x))
            used += tm0_list + tm1_list + [makeMate(x) for x in tm0_list+tm1_list]
    if is_closer and signal < 0:
        return True
    elif is_closer and signal > 0:
        return False
    elif not(is_closer) and signal < 0:
        return False
    elif not(is_closer) and signal > 0:
        return True 
    else:
        return ValueError("signal must be non-zero")




class MarkerMiniClass:
    def __init__(self,G,gene_assoc,m_id,car_id,locus_range,current_nb):
        self.reference_graph = G
        self.nb_candidate = current_nb
        self.gene_assoc = gene_assoc
        self.m_id = m_id
        self.car_id = car_id
        self.locus_range = locus_range
    def __lt__(self,m):
        if self.locus_range != m.locus_range:
            return self.locus_range < m.locus_range
        return nb_search(self.reference_graph,self.gene_assoc,self.m_id,m.m_id,self.car_id,self.nb_candidate,-1)
    def __gt__(self,m):
        if self.locus_range != m.locus_range:
            return self.locus_range > m.locus_range
        return nb_search(self.reference_graph,self.gene_assoc,self.m_id,m.m_id,self.car_id,self.nb_candidate,1)
    def __le__(self,m):
        return self.locus_range <= m.locus_range
    def __ge__(self,m):
        return self.locus_range >= m.locus_range
    def __eq__(self,m):
        return self.locus_range == m.locus_range
    def __ne__(self,m):
        return self.locus_range != m.locus_range

def createHeap(G,nb_list,gene_assocs,car_id,current_nb):
    new_list = []
    for v in nb_list:
        new_list.append(MarkerMiniClass(G,gene_assocs,v,car_id,gene_assocs[v][car_id],current_nb))
    return BinaryHeap(new_list)

def chooseNextNeighbour(nb_heap,profiles):
    v_item = nb_heap.pop()
    while profiles[v_item.m_id] <= 0:
#        if v_item.m_id == ('171',0):
#            print v_item.m_id
        v_item = nb_heap.pop()
    return v_item.m_id

def contractWindow(u,gene_assocs,profiles,car_id):
    contractable = False        
    if profiles[u] == 0:
        gene_assocs[u][car_id] = [sys.maxint]
    elif profiles[u] == 1:
        for i in range(len(gene_assocs[u][car_id])-1):
            if abs(gene_assocs[u][car_id][i+1]-gene_assocs[u][car_id][i]) > 1:
                contractable = True
                break
        if contractable:
            gene_assocs[u][car_id] = gene_assocs[u][car_id][i+1:]
        else:
            k = len(gene_assocs[u][car_id])>>1
            if k == 0:
                gene_assocs[u][car_id] = [sys.maxint]
            else:
                gene_assocs[u][car_id] = gene_assocs[u][car_id][k:]
                        
def addToPath(G,path_list,u,profiles,gene_assocs,car_id):
    path_list.append(u)
    path_list.append(makeMate(u))
    profiles[u] -= 1
    profiles[makeMate(u)] -= 1
    candidate_nbs = [v for v in G[makeMate(u)] \
            for k in G[makeMate(u)][v] \
            if (G[makeMate(u)][v][k]['color'] != 'b' and \
            profiles[v] > 0)]
    contractWindow(u,gene_assocs,profiles,car_id)
    contractWindow(makeMate(u),gene_assocs,profiles,car_id)
    return candidate_nbs

def redLen(G,v):
    return len(redNbhood(G,v))

def redNbhood(G,v):
    return [u for u in G[v] for k in G[v][u] if G[v][u][k]['color'] == 'r']

def findMarkerOrder(G,mults,gene_assocs,ends,car_id):
    profiles = copy.deepcopy(mults)
    used = {}
    for e in ends:
        used[e] = False
    paths = []
    while any([not(used[u]) for u in used]):
        for e in ends:
            if not(used[e]):
                end0 = e
                break
        current_path = []            
        current_path.append(end0)
        current_path.append(makeMate(end0))
        used[end0] = True
        profiles[end0] -= 1
        profiles[makeMate(end0)] -= 1
        nbhood = [v for v in G[makeMate(end0)] \
                for k in G[makeMate(end0)][v] \
                if G[makeMate(end0)][v][k]['color'] != 'b']
        current_vx = current_path[-1]
        while nbhood != [] and any([profiles[v] > 0 for v in nbhood]):
            if len(nbhood) > 1 and \
                    makeMate(current_vx) in nbhood and \
                    all([profiles[x] == 0 for x in G[makeMate(current_vx)]]):
                profiles[current_vx] -= 1
                profiles[makeMate(current_vx)] -= 1
            elif len(nbhood) == 1:
                current_vx = nbhood.pop()
                nbhood = addToPath(G,current_path,current_vx,profiles,gene_assocs,car_id)
            elif len(nbhood) > 1:
                to_print = False
                if current_vx == ('95',1) or current_vx == ('96',1):
                    to_print = False 
                if to_print:
                    for x in nbhood:
                        print x,profiles[x],gene_assocs[x][car_id]
                if len(nbhood) == 2 and current_vx in nbhood and \
                        makeMate(current_vx) in nbhood and \
                        mults[current_vx] > profiles[current_vx]:
                    current_vx = makeMate(current_vx)
                else:                
                    current_heap = createHeap(G,nbhood,gene_assocs,car_id,current_path[-1]) 
                    current_vx = chooseNextNeighbour(current_heap,profiles)
                if to_print:
                    print "Chosen",current_vx
                nbhood = addToPath(G,current_path,current_vx,profiles,gene_assocs,car_id)
                if to_print:
                    print "New choose set",nbhood
            if current_path[-1] in used:
                used[current_path[-1]] = True
            else:
                if max(gene_assocs[current_path[0]][car_id]) > \
                        max(gene_assocs[current_path[-1]][car_id]):
                    current_path = current_path[::-1]
        if current_path != [] and paths != []:
            paths.append('|')
        paths += current_path                
    return paths

def findAltCycle(G,v,traversal=[],positions={}):
    nbhood = G[v]
    future_traversals = []
    for u in G[v]:
        if u in positions:
            return traversal,positions[u],positions[v]
        new_positions = copy.deepcopy(positions)
        new_traverse = traversal[:]
        new_traverse.append(v)
        new_traverse.append(makeMate(v))
        new_positions
        new_traverse = findAltCycle() 
    stack = [v]
    prev_branching = 0
    while stack != [v]:
        current_vertex = stack.pop()
        current_mate = makeMate(current_vertex)
        traversed.append(current_vertex,current_mate)
        if len(G[current_mate]) == 1:
            traversed = traversed[:prev_branching]
        elif len(G[current_mate]) > 2:
            prev_branching = len(traversed) 
        for u in G[current_mate]:
            if u != current_vertex and u not in traversed:
                if u == makeMate(v):
                    return traversed
                elif u not in traversed:
                    stack.append(u)
    return []


def areAllDiffWindows(window_limits):
    is_diff = True
    for i in window_limits:
        for j in window_limits:
            if i != j:
                if window_limits[i] == window_limits[j]:
                    return False
    return True                

def breadthFirst(G,start_vertex,marker_windows):
    search_array = [start_vertex]
    populate_array = [v for v in G[start_vertex] if v != makeMate(start_vertex)]
    branches = [(start_vertex,v,makeMate(v),) for v in G[start_vertex] if v != makeMate(start_vertex)]
    limiting_windows = {}
    for i in range(len(branches)):
        limiting_windows[i] = \
                [max([x for x in marker_windows[v] for v in branches[i]]),min([x for x in marker_windows[v] for v in branches[i]])]                
    extendable = True                
    while not(areAllDiffWindows) and extendable:
        add_to_branches = {}                
        extend_array = []
        for i in range(len(branches)):
            nbhood = G[branches[i][-1]]
            if len(nbhood) == 0:
                extend_array.append(False)
                continue
            extend_array.append(True)
            add_to_branches[i] = []
            if nbhood[0] not in branches[i] \
                and makeMate(nbhood[0]) not in branches[i]:
                branches[i] += (nbhood[0],makeMate(nbhood[0]))
                limiting_windows[i] = \
                        [max([x for x in marker_windows[v] for v in branches[i]]),min([x for x in marker_windows[v] for v in branches[i]])]
            if len(nbhood) > 2:
                for j in range(1,len(nbhood)):
                    if nbhood[j] not in branches[i] and makeMate(nbhood[j]) not in branches[i]:
                        add_to_branches.append(branches[i]+(nbhood[j],makeMate(nbhood[j])))
        for i in add_to_branches:
            for j in range(len(add_to_branches[i])):
                b = add_to_branches[i][j]
                branches.append(b)
                limiting_windows[len(branches)-1] = \
                        [max([x for x in marker_windows[v] for v in branches[len(branches)-1]]),\
                        min([x for x in marker_windows[v] for v in branches[len(branches)-1]])]
        extendable = any(extend_array)    
    return branches,limiting_windows



def componentTraverser(G_comp,car_id,gene_assocs,mults):
    win_ends = {}
    for g in G_comp:
        if g not in win_ends:
            win_ends[g] = None
        g_ends = findEnds(g,gene_assocs,mults,car_id)
        win_ends[g] = g_ends
    sorted_components = sorted(G_comp,key=lambda x: max([min(gene_assocs[v][car_id]) for v in win_ends[x]]))
    join_components = []
    for i in range(len(sorted_components)):
        g = sorted_components[i]
        join_components += findMarkerOrder(g,mults,gene_assocs,win_ends[g],car_id)
        if i < len(sorted_components)-1:
            join_components.append('|')
    return join_components
