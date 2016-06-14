#!/usr/bin/env

import re
import sys
import copy

# Class for reading and storing NHX tree
class NHXTree:
    node_info = ''
    node_name = ''
    branch_length = 0
    children_list = []
    parent = None 
    def __init__(self, 
                 node_info,
                 children_list,length = None):
        self.node_info = node_info
        self.children_list = children_list
        self.node_name = re.search('^([^:\[;]*)[:\[;]?.*$',node_info).group(1)
        if length == None:
            pos = self.node_info.find(':')+1
            k = pos
            while k < len(self.node_info) and (self.node_info[k].isdigit() or self.node_info[k] == '.'):
                k += 1
            if node_info[pos:k] != '':    
                self.branch_length = float(node_info[pos:k])
        else:
            self.branch_length = length
    def get_children(self):
        return copy.copy(self.children_list)
    def __str__(self):
        nhx_string = ''
        if len(self.children_list) > 0:
            nhx_string = '('
            for child in self.children_list:
                nhx_string += str(child)
                nhx_string += ','
            nhx_string = nhx_string[:-1]
            nhx_string += ')'
        if self.node_info == '':
            nhx_string += ';'
            return nhx_string
        nhx_string += self.node_info
        return nhx_string
    def set_parent_root(self):
        for i in range(len(self.children_list)):
            self.children_list[i].set_parent(self)
    def set_parent(self,p):
        self.parent = p
        for i in range(len(self.children_list)):
            self.children_list[i].set_parent(self)
    def rename(self,label='',counter=0):
        if len(self.children_list) > 0:
            if self.node_info[0]==':' and label!='':
                self.node_info=label+'_'+str(counter)+self.node_info
                counter+=1
            for child in self.children_list:
                counter=child.rename(label,counter)
        return counter

    def search(self, name):
        """Returns a list of subtrees with the given root."""
        if self.node_name == name:
            return [self]
        results = []
        for child in self.children_list:
            subtree = child.search(name)
            if subtree:
                results += subtree
        return results

    def leaves(self):
        """Returns list of leaves"""
        if len(self.children_list)==0:
            return [self]
        list_of_leaves=[]
        for child in self.children_list:
            list_of_leaves += child.leaves()
        return list_of_leaves
    
    def dollo_parsimony_items(self,leaf_items):
        internal_labels = self.bottom_up_items(leaf_items)
        return internal_labels
    
    def bottom_up_dollo(self,leaf_labels,to_print=False):
        all_labels = {}
        if len(self.children_list) == 0 and self.node_name in leaf_labels:
            all_labels = {self:leaf_labels[self.node_name]}
            return all_labels
        elif len(self.children_list) == 0 and self.node_name not in leaf_labels:
            all_labels = {self:'0'}
            return all_labels
        candidate = '2'
        count = 0
        uncertain_count = 0
        for child in self.children_list:
            all_labels.update(child.bottom_up_dollo(leaf_labels))
            if all_labels[child] == '1' or all_labels[child] == '2':
                count += 1
        if count > 1:# len(self.children_list):
            candidate = '1'
        elif count == 0:
            candidate = '0'
        else:
            candidate = '2'
        all_labels[self] = candidate
        return all_labels

    def top_down_dollo(self,new_labels,to_set = '0'):
        if self.parent == None:
            count = 0
            for c in self.children_list:
                if new_labels[c] == '1' or new_labels[c] == '2':
                    count += 1
            if count > 1:
                to_set = '1'
            else:
                to_set = '0'
        if new_labels[self] == '2' and to_set =='1':
            new_labels[self] = '1'
        elif new_labels[self] == '2' and to_set =='0':
            new_labels[self] = '0'
        for child in self.children_list:
            child.top_down_dollo(new_labels,new_labels[self])
            
    def bottom_up_items(self,leaf_labels):
        # If leaf, return leaf label
        if len(self.children_list) == 0 and self.node_name in leaf_labels:
            all_labels = {self:leaf_labels[self.node_name]}
            return all_labels
        elif len(self.children_list) == 0 and self.node_name not in leaf_labels:
            all_labels = {self:0}
            return all_labels
        else:
            current_max = -1
            all_labels = {}
            for child in self.children_list:
                all_labels.update(child.bottom_up_items(leaf_labels))
                if all_labels[child] > current_max:
                    current_max = all_labels[child]
            if current_max > 0:
                all_labels[self] = current_max
            else:
                all_labels[self] = 0
        return all_labels

    def fitch_parsimony(self,leaf_labels,choice):
        internal_labels = self.bottom_up_fitch(leaf_labels)
        internal_labels = self.top_down_fitch(internal_labels,choice)
        return internal_labels

    # Bottom-up phase
    def bottom_up_fitch(self,leaf_labels):
        # If leaf, return leaf label
        if len(self.children_list) == 0 and self.node_name in leaf_labels:
            all_labels = {self:leaf_labels[self.node_name]}
            return all_labels
        # If leaf, but unassigned, return error
        if len(self.children_list) == 0 and self.node_name not in leaf_labels:
            all_labels = {self:'0/1'}
            return all_labels
        # Method for internal nodes
        else:
            # Default assignment
            combine = '0/1'
            all_labels = {}
            # Fitch parsimony bottom-up phase: get assignments from children
            for child in self.children_list:
                all_labels.update(child.bottom_up_fitch(leaf_labels))
                if combine == '0' and all_labels[child] == '0':
                    combine = '0'
                elif combine == '1' and all_labels[child] == '1':
                    combine = '1'
                elif combine == '0' and all_labels[child] == '1':
                    combine = '0/1'
                elif combine == '1' and all_labels[child] == '0':
                    combine = '0/1'
                elif combine == '0/1' and all_labels[child] == '1':
                    combine = '1'            
                elif combine == '0/1' and all_labels[child] == '0':
                    combine = '0'
                elif combine == '0/1' and all_labels[child] == '0/1':
                    combine = '0/1'
                elif combine == '0' and all_labels[child] == '0/1':
                    combine = '0'            
                elif combine == '1' and all_labels[child] == '0/1':
                    combine = '1'
            all_labels[self] = combine
            # If the node has already been assigned 0 or 1 in the initial assignment, revert.
            if self.node_name in leaf_labels and leaf_labels[self.node_name] not in ['0/1','2']:
                all_labels[self] = leaf_labels[self.node_name]
            return all_labels


    # Top-down phase    
    def top_down_fitch(self,labels,choice):
        # Make choice for 0/1
        if labels[self] == '0/1':
            labels[self] = choice
        # If leaf, return labels
        if len(self.children_list) == 0:
            return labels
        # Propogate labels to children
        for child in self.children_list:
            if labels[self] == '0':
                    if labels[child] == '0/1':
                        labels[child] = '0'
                    child.top_down_fitch(labels,0)
            if labels[self] == '1':
                for child in self.children_list:
                    if labels[child] == '0/1':
                        labels[child] = '1'
                    child.top_down_fitch(labels,1)
        return labels
    
    def dollo_parsimony(self,leaf_labels,choice,to_print=False):
        internal_labels = self.bottom_up_dollo(leaf_labels,to_print)
        self.top_down_dollo(internal_labels)
        return internal_labels
    
    def twice_conserved(self,leaf_labels,dummy = '1'):
        count = 0
        for l in leaf_labels:
            if leaf_labels[l] == '1':
                count += 1
        if count > 1:
            return '1'
        else:
            return '0'

def sankoff_rousseau(T,leaf_labels,states,printer=False):
    orig_labels = initializer(T,leaf_labels,states)
    new_labels = sankoff_rousseau_bottom_up(T,orig_labels,states)
    best_labels = sankoff_rousseau_top_down(T,new_labels,-1,printer)
    return best_labels


def initializer(T,leaf_labels,states):
    init_dict = {}
    for l in T.leaves():
        init_dict[l] = {}
        if l.node_name in leaf_labels:
            init_dict[l][leaf_labels[l.node_name]] = 0
            for s in states:
                if s != leaf_labels[l.node_name]:
                    init_dict[l][s] = sys.maxint
        else:
            init_dict[l][0] = 0
            for s in states:
                if s != 0:
                    init_dict[l][s] = sys.maxint
    return init_dict


def sankoff_rousseau_bottom_up(T,orig_labels,states,printer = False):
    current_labels = {}
    if len(T.children_list) == 0:# orig_labels:
        current_labels[T] = orig_labels[T]#.update({T:orig_labels[T]})
        return current_labels
    else:
        for c in T.children_list:
            current_labels.update(sankoff_rousseau_bottom_up(c,orig_labels,states,printer))
        current_labels[T] = {}     
        for s in states:
            total = 0
            for c in T.children_list:
                minimizer = s
                minimum_score = sys.maxint
                for t in states:
                    if minimum_score > current_labels[c][t] + abs(s-t):
                        minimizer = t
                        minimum_score = current_labels[c][t] + abs(s-t)
                total += minimum_score
            current_labels[T][s] = total                
        return current_labels

def sankoff_rousseau_top_down(T,internal_labels,prev_label=-1,printer=False):
    factor = 1
    if prev_label < 0:
        factor = 0
    best_labels = {}        
    best_score = sys.maxint
    local_best = None
    local_dict = internal_labels[T]
    if printer:
        print T,local_dict
    for t in local_dict:
        if local_dict[t]+(factor*abs(prev_label-t)) < best_score or \
            (local_dict[t]+(factor*abs(prev_label-t)) == best_score and \
            t > local_best):
            local_best = t
            best_score = local_dict[t]+(factor*abs(prev_label-t))
    best_labels[T] = local_best
    for child in T.children_list:
        best_labels.update(sankoff_rousseau_top_down(child,internal_labels,local_best,printer))
    return best_labels


def rerooted(tree,leaf):
    l = leaf.branch_length
    T = NHXTree('',[NHXTree(leaf.node_info,[],0),makeRoot(tree,leaf.parent,leaf,l)])
    return T

def makeRoot(tree,current_node,to_remove,length):
    # Name of node to add
    # Name of node to make child, i.e. parent
    make_child = current_node.parent
    # Next branch length for new child
    new_length = current_node.branch_length
    # New child list of current node 
    new_children_list = current_node.get_children()
    # Removing previously added node, which is now a parent
    new_children_list.remove(to_remove)
    # If the next thingie is a dud
    if make_child == None:
        return NHXTree(current_node.node_info,new_children_list,length) 
    else:
        new_children_list.append(makeRoot(tree,make_child,current_node,new_length))
        T = NHXTree(current_node.node_info,new_children_list,length)
        return T 

def get_children_strings(string):
    """Returns a list of the strings for each child."""
    par_count = 0
    children_list = []
    start = 0

    for i in range(len(string)):
        if string[i] == '(':
            par_count += 1
        elif string[i] == ')':
            par_count -= 1
        elif string[i] == ',' and par_count == 0:
            children_list.append(string[start:i])
            start = i+1
        elif i == len(string)-1:
            children_list.append(string[start:])
    return children_list

def parse_nhx(string,parent=None):
    """Parses the NHX tree.

    Returns a NHXTree data structure.
    
    """
    end = len(string)-1
    # Get the root
    while end > -1 and string[end] != ')':
        end -= 1
    root_info = string[end+1:]
    # Parse the children
    children_list = []
    if end > 1:
        children_strings = get_children_strings(string[1:end])
        if len(children_strings) > 0:
            for child in children_strings:
                children_list.append(parse_nhx(child))

    # Build the data structure
    root = NHXTree(root_info,children_list)
    return root



def queryLCA(tree,leaf_set):
    call_stack = [tree]
    lca = None
    while call_stack != []:
        current_ancestor = call_stack.pop()
        bool_ancestor = True
        for i in range(len(leaf_set)):
            bool_ancestor = bool_ancestor and isAncestral(current_ancestor,leaf_set[i])
        if bool_ancestor:
            lca = current_ancestor
            call_stack = call_stack + current_ancestor.children_list
    return lca

def isAncestral(tree,leaf):
    tree.set_parent_root()
    current_parent = leaf
    while current_parent != None:
        current_parent = current_parent.parent
        if current_parent == tree:
            return True
    return False


def prune_tree(tree,label_list):
    if tree.name not in label_list and tree.parent != None:
        p = tree.parent
        b = tree.branch_length
        c_list = tree.chilren_list
        p.children_list.remove(tree)
        del tree
        for child in c_list:
            child.branch_length += b
            p.children_list.append(child)

    
        

