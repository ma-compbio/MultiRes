#!/usr/bin/env
# Module to facilitate manipulating adjacencies
from marker_string import *
import random
import kmp
import lcs

def tupleToDict(tuple_list,add_mates=True):
    hash_table = {}
    for t in tuple_list:
        t0,t1 = t
        if t0 not in hash_table:
            hash_table[t0] = []
            if add_mates:
                hash_table[t0].append(makeMate(t0))
        if makeMate(t0) not in hash_table:            
            hash_table[makeMate(t0)] = []
            if add_mates:
                hash_table[makeMate(t0)].append(t0)
        if t1 not in hash_table:
            hash_table[t1] = []
            if add_mates:
                hash_table[t1].append(makeMate(t1))
        if makeMate(t1) not in hash_table:            
            hash_table[makeMate(t1)] = []
            if add_mates:
                hash_table[makeMate(t1)].append(t1)
        if t1 not in hash_table[t0]:
            hash_table[t0].append(t1)
        if t0 not in hash_table[t1]:
            hash_table[t1].append(t0)
    return hash_table                        

def dfs(adj_dict):
    visited = set([])
    components = []
    while len(visited) != len(adj_dict):
        stack = [random.choice(adj_dict.keys())]
        while stack[-1] in visited:
            stack = [random.choice(adj_dict.keys())]
        current_component = [stack[-1]]
        visited.add(stack[-1])
        while len(stack) > 0:
            current_vertex = stack.pop()
            for v in adj_dict[current_vertex]:
                if v not in visited:
                    stack.append(v)
                    visited.add(v)
                    current_component.append(v)
        components.append(current_component)
    return components        

def findConnectedComponents(adjacencies):
    adj_dict = tupleToDict(adjacencies)
    connected_components = dfs(adj_dict)
    #print "CC:",len(connected_components)
    #print connected_components
    return connected_components

