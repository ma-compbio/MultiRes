#!/usr/bin/env
import sys

# Naive suffix array construction.
# SL classification of characters. True for S-type, False for L-type.

def classifySL(string):
    sl_array = [False for x in string]
    lms_dict = []
    end = len(string)
    for i in range(len(string)-1,-1,-1):
        if i == len(string)-1:
            sl_array[i] = False
        else:
            if(string[i] < string[i+1]) or (string[i] == string[i+1] and sl_array[i+1] == True):
                sl_array[i] = True
            else:
                sl_array[i] = False
                if sl_array[i+1] == True:
                    lms_dict.append((i+1,end))
                    end = i+2
    lms_dict.reverse()                        
    return sl_array,lms_dict

def compareLMSsubstrings(string,pos1,pos2,sl_array):
    s1,t1 = pos1[0]
    s2,t2 = pos2[0]
    while s1 < t1 and s2 < t2:
        if string[s1] > string[s2]:
            return pos2
        elif string[s1] < string[s2]:
            return pos1
        elif string[s1] == string[s2]:
            if sl_array[s1] == True and sl_array[s2] == False:
                return pos1
            elif sl_array[s2] == True and sl_array[s1] == False:
                return pos2
            else:
                s1 += 1
                s2 += 1
    return None

def suffixArray(string):
    temp_array = []
    for i in range(len(string)):
        temp_array.append((string[i:],i))
    temp_array.sort(key=lambda x: x[0])
    return [temp_array[i][-1] for i in range(len(temp_array))]

def getInverse(suffix_array):
    inv = [0 for i in range(len(suffix_array))]
    for i in range(len(suffix_array)):
        inv[suffix_array[i]] = i
    return inv

def lcpArray(string,suffix_array):
    inv_array = getInverse(suffix_array)
    lcp_array = [-1 for i in range(len(suffix_array))]
    for i in range(len(string)):
        l = 0
        k = inv_array[i]
        if k == 0:
            continue
        j = suffix_array[k-1]
        while i+l < len(string) and j+l < len(string) and string[i+l] == string[j+l]:
            l += 1
        lcp_array[k] = l
    return lcp_array

def longestSubstring(string1,string2,end_points = False):
    new_string = string1+string2
    part_index = len(string1)
    suff_array = suffixArray(new_string)
    inv_array = getInverse(suff_array)
    lcp_array = lcpArray(new_string,suff_array)
    score = -1
    pos1 = len(string1+string2)
    pos0 = -1
    for i in range(0,len(new_string)):
        suff_index = inv_array[i]
        other_index = suff_array[suff_index-1]
        if min(i,other_index) >= part_index or \
                max(i,other_index) < part_index:
            continue                    
        new_score = lcp_array[suff_index]
        if min(i,other_index) < part_index and \
                min(i,other_index)+new_score >= part_index:
                new_score = part_index-min(i,other_index)
        if new_score > 0 and new_score >= score:
            if end_points and \
                    (pos0 < min(i,other_index)):
                pos0 = min(i,other_index)
            else:
                pos0 = min(i,other_index)
            if end_points and pos1 > max(i,other_index):
                pos1 = max(i,other_index)
            else:
                pos1 = max(i,other_index)
            score = new_score
    return score,pos0,pos1-part_index

# Edit-distance and and optimal alignment:
#   default scores = (1,1,1)
def editDistance(string1,string2,mismatch=1,insert=1,delete=1):
    dp_array = [[0 for j in range(len(string2)+1)] for i in range(len(string1)+1)]
    for i in range(1,len(string1)+1):
        dp_array[i][0] = delete*i
    for j in range(1,len(string2)+1):
        dp_array[0][j] = insert*j
    for i in range(1,len(string1)+1):
        for j in range(1,len(string2)+1):
            dp_array[i][j] = min(insert+dp_array[i][j-1],delete+dp_array[i-1][j],\
                    matchScore(string1[i-1],string2[j-1],mismatch)+dp_array[i-1][j-1])
    return editBacktrackIter(dp_array,string1,string2,mismatch,insert,delete),dp_array[-1][-1]

def localEditDistance(string1,string2,mismatch=1,insert=1,delete=1):
    dp_array = [[0 for j in range(len(string2)+1)] for i in range(len(string1)+1)]
    for i in range(1,len(string1)+1):
        dp_array[i][0] = 0 
    for j in range(1,len(string2)):
        dp_array[0][j] = 0
    for i in range(1,len(string1)+1):
        for j in range(1,len(string2)+1):
            dp_array[i][j] = min(insert+dp_array[i][j-1],delete+dp_array[i-1][j],\
                    matchScore(string1[i-1],string2[j-1],mismatch)+dp_array[i-1][j-1])
    cost = len(string1)+len(string2)+1
    start_local = -1
    for i in range(len(string1),0,-1):
        if dp_array[i][len(string2)] < cost:
            cost = dp_array[i][len(string2)]
            start_local = i
    return editBacktrackIter(dp_array,string1,string2,mismatch,insert,delete,start_local),cost


def matchScore(element1,element2,cost):
    if element1 == element2:
        return 0
    else:
        return cost

def editBacktrack(dp_array,string1,string2,i,j,mismatch,insert,delete):
    if i == 0 and j == 0:
        return [(string1[i],string2[j])]
    elif i == 0 or j == 0:
        to_return = [(string1[i],string2[j])]
        if i > j:
            for k in range(i-1,-1,-1):
                to_return.append((string1[k],'_'))
        else:
            for k in range(j,-1,-1):
                to_return.append(('_',string2[k]))
        to_return.reverse()
        return to_return
    if dp_array[i][j] == dp_array[i-1][j-1]+matchScore(string1[i],string2[j],mismatch):
        return editBacktrack(dp_array,string1,string2,i-1,j-1,mismatch,insert,delete)+[(string1[i],string2[j])]
    elif dp_array[i][j] == dp_array[i-1][j] + delete:
        return editBacktrack(dp_array,string1,string2,i-1,j,mismatch,insert,delete)+[(string1[i],'_')]
    elif dp_array[i][j] == dp_array[i][j-1] + insert:
        return editBacktrack(dp_array,string1,string2,i,j-1,mismatch,insert,delete)+[('_',string2[j])]

def editBacktrackIter(dp_array,string1,string2,mismatch,insert,delete,start_local=-1):
    i = len(string1)
    j = len(string2)
    if start_local != -1:
        i = start_local
    alignment = []
    while i > 0 or j > 0:
        if i > 0  and j > 0 and \
                dp_array[i][j] == dp_array[i-1][j-1]+matchScore(string1[i-1],string2[j-1],mismatch):
            alignment.append((string1[i-1],string2[j-1]))
            i -= 1
            j -= 1
        elif dp_array[i][j] == dp_array[i-1][j] + delete:
            alignment.append((string1[i-1],'_'))
            i -= 1
        elif dp_array[i][j] == dp_array[i][j-1] + insert:
            alignment.append(('_',string2[j-1]))
            j -= 1
        elif i == 0 or j == 0:
            if i == j:
                i -= 1
                j -= 1
            elif i == 0 and j > 0:
                alignment.append(('_',string2[j-1]))
                j -= 1
            elif j == 0 and i > 0:
                alignment.append((string1[i-1],'_'))
                i -= 1
    alignment = alignment[::-1]                
    if start_local != -1:
        if start_local+1 < len(string1)-1:
            for k in range(start_local+1,len(string1)):
                alignment.append((string1[k],'_'))
    return alignment

def isSubsequence(string0,string1):
    if len(string0) < len(string1):
        return -1 
    k = 0
    score = 0
    for i in range(len(string1)):
        while k < len(string0) and\
                string0[k] != string1[i]:
                    k += 1
        if k < len(string0):
            score += 1
        k += 1
    return score 
def contains_sublist(lst, sublst):
    n = len(sublst)
    return any((sublst == lst[i:i+n]) for i in range(len(lst)-n+1))



# Longest common subsequence for 2 strings

def dpLongestCommonSubsequence(string1,string2):
    dp_array = [[0 for j in range(len(string2))] for i in range(len(string1))]
    for i in range(len(string1)):
        for j in range(len(string2)):
            if string1[i] != string2[j]:
                dp_array[i][j] = max(dp_array[i][j-1],dp_array[i][j-1])
            else:
                dp_array[i][j] += 1+dp_array[i-1][j-1]
    subsequence = ""# subsequence_backtrack(dp_array,string1,string2,len(string1)-1,len(string2)-1)
    return subsequence,dp_array[-1][-1]

def subsequence_backtrack(dp_array,string1,string2,i,j):
    if i == -1 or j == -1:
        return ""
    elif string1[i] == string2[j]:
        return subsequence_backtrack(dp_array,string1,string2,i-1,j-1)+string1[i]
    elif dp_array[i][j-1] > dp_array[i-1][j]:
        return subsequence_backtrack(dp_array[i][j-1],string1,string2,i,j-1)
    else:
        return subsequence_backtrack(dp_array[i][j-1],string1,string2,i-1,j)

# LCS algorithm implmented as a dynamic program; 
# this is inefficient, compared to using suffix trees/arrays.
# But as long as strings are small, it should suffice.

def dpLongestCommonSubstring(string1,string2):
    dp_array = [[0 for j in range(len(string2))] for i in range(len(string1))]
    max_score = 0
    big_string_indices = []
    for i in range(len(string1)):
        for j in range(len(string2)):
            if string1[i] == string2[j]:
                if i == 0 or j == 0:
                    dp_array[i][j] = 1
                else:
                    dp_array[i][j] = 1 + dp_array[i-1][j-1]
                if dp_array[i][j] > max_score:
                    max_score = dp_array[i][j]
                    big_string_indices = [i-max_score+1,i+1]
            else:
                dp_array[i][j] = 0
    return max_score,big_string_indices


def reverser(string):
    reversed_string = []
    for x in reversed(string):
        if x.isspace() or len(x) == 0:
            continue
        if x[0] == '+':
            reversed_string.append('-'+x[1:])
        elif x[0] == '-':
            reversed_string.append('+'+x[1:])
    return reversed_string

def reducer(t0,t1):
    reduced = []
    if '_h' in t0:
        reduced.append('+'+t0[:-2])
    elif '_t' in t0:
        reduced.append('-'+t0[:-2])
    if '_h' in t1:
        reduced.append('-'+t1[:-2])
    elif '_t' in t1:
        reduced.append('+'+t1[:-2])
    return reduced

def expander(t,is_num=True):
    if t[0] != '-':
        if is_num:
            id =[2*int(t)-1,2*int(t)]
            return ' '.join([str(x) for x in id])
        else:
            return t[1:]+'_t'+' '+t[1:]+'_h'
    else:
        if is_num:
            id =[2*abs(int(t)),2*abs(int(t))-1]
            return ' '.join([str(x) for x in id])
        else:
            return t[1:]+'_h'+' '+t[1:]+'_t'

def bigExpander(interval,is_num=True):
    fill_in = []
    for x in interval:
        fill_in.append(expander(x,is_num))
    return fill_in# ' '.join(fill_in)

def parseContainments(blocks,in_species,out_species=[],adjacency_weights={}):
    parse_containments = {}
    span_id = 0
    convolution_factor = 0.9
    for b in blocks:
        counts,associations,out_associations = getCounts(blocks[b],in_species,out_species,to_print)
        do_other = False 
        if not(do_other):
            list_containments,weights,sps = sublist_checker(blocks[b],in_species,out_species)
            if len(adjacency_weights) == 0:
                for i in range(len(list_containments)):
                    print str(span_id)+'|'+str(weights[i])+';'+','.join(sps[i])+':'+' '.join(bigExpander(list_containments[i]))
                    span_id += 1
            else:
                content = tuple(b.split(' '))
                contrib = 0.
                if content in adjacency_weights:
                    contrib = adjacency_weights[content]
                elif (content[-1],content[0]) in adjacency_weights:
                    contrib = adjacency_weights[(content[-1],content[0])]
                for i in range(len(list_containments)):
                    weight = weights[i]
                    if contrib != 0:
                        weight = convolution_factor*(contrib) + (1-convolution_factor)*weight
                    print str(span_id)+'|'+str(weight)+';'+','.join(sps[i])+':'+' '.join(bigExpander(list_containments[i]))
                    span_id += 1
            continue
        sorted_counts = sorted([x for x in counts],key=lambda x: counts[x],reverse=True)
        for c in sorted_counts:
            if counts[c] >= (len(in_species)>>1) and counts[c] < len(in_species):
                if len(out_species) == 0:
                    weight = 0.
                    if len(adjacency_weights) == 0:
                        weight = 1.*counts[c]/(len(in_species))
                    else:
                        content = tuple(b.split(' '))
                        contrib = 0.
                        weight = 1.*counts[c]/(len(in_species))
                        if content in adjacency_weights:
                            contrib = adjacency_weights[content]
                        elif (content[-1],content[0]) in adjacency_weights:
                            contrib = adjacency_weights[(content[-1],content[0])]
                        if contrib != 0:
                            weight = convolution_factor*(contrib) + (1-convolution_factor)*weight
#                        else:
#                            weight = 1.*counts[c]/(len(in_species))
                    if len(parse_containments[b]) > 2:
                        print str(span_id)+'|'+str(weight)+';'+c+':'+' '.join(bigExpander(parse_containments[b]))
                        span_id += 1
                else:        
                    support = 0
                    for i in range(len(associations[c])-1):
                        found = False
                        for other_sp in out_associations:
                           for k in range(len(out_associations[other_sp])-1):
                               if out_associations[other_sp][k] == associations[c][i] and \
                                       out_associations[other_sp][k+1] == associations[c][i+1]:
                                           support +=1
                                           found = True
                                           break
                           if found:
                               break
                    if support >= ((len(associations[c])-1)>>1):
                        if len(associations[c]) <= 1:
                            break
                        parse_containments[b] = associations[c]
                        weight = 0.
                        if len(adjacency_weights) == 0:
                            weight = 1.*counts[c]/(len(in_species))
                        else:
                            content = tuple(b.split(' '))
                            contrib = 0.
                            weight = 1.*counts[c]/(len(in_species))
                            if content in adjacency_weights:
                                contrib = adjacency_weights[content]
                            elif (content[-1],content[0]) in adjacency_weights:
                                contrib = adjacency_weights[(content[-1],content[0])]
                            if contrib != 0:
                                weight = convolution_factor*(contrib) + (1-convolution_factor)*weight
                        if len(parse_containments[b]) > 2:
                            print str(span_id)+'|'+str(weight)+';'+c+':'+' '.join(bigExpander(parse_containments[b]))
                            span_id += 1
                            break
            elif counts[c] == len(in_species):
                if len(associations[c]) <= 1:
                    break
                parse_containments[b] = associations[c]
                weight = 1.
                if len(adjacency_weights) == 0:
                    weight = 1.*counts[c]/(len(in_species))
                else:
                    content = tuple(b.split(' '))
                    contrib = 0.
                    if content in adjacency_weights:
                        contrib = adjacency_weights[content]
                    elif (content[-1],content[0]) in adjacency_weights:
                        contrib = adjacency_weights[(content[-1],content[0])]
                    if contrib != 0:
                        weight = convolution_factor*(contrib) + (1-convolution_factor)*weight
                if len(parse_containments[b]) > 2:
                    print str(span_id)+'|'+str(weight)+';'+c+':'+' '.join(bigExpander(parse_containments[b]))
                    span_id += 1
                    break_out = True
                    break


# Check largest sublist

def sublist_checker(block,in_species,out_species):
    lists_contained = []
    weights = []
    species = []
    for b1 in block:
        sp1 = b1[:b1.find('.')]
        if sp1 in in_species:
            for b2 in block:
                sp2 = b2[:b2.find('.')]
                if b1 != b2 and sp2 in out_species:
                    if contains_sublist(block[b1],block[b2]) and block[b2] not in lists_contained:
                        if len(block[b2]) > 2:
                            lists_contained.append(block[b2])
                            weights.append(2)
                            species.append([sp1,sp2])
                    elif contains_sublist(block[b2],block[b1]) and block[b1] not in lists_contained:
                        if len(block[b1]) > 2:
                            lists_contained.append(block[b1])
                            weights.append(2)
                            species.append([sp1,sp2])
                    elif block[b1] in lists_contained:
                        for i in range(len(lists_contained)):
                            if lists_contained[i] == block[b1]:
                                if sp1 not in species[i]:
                                    weights[i] += 1
                                    species[i].append(sp1)
                                break
                    elif  block[b2] in lists_contained:
                        for i in range(len(lists_contained)):
                            if lists_contained[i] == block[b2]:
                                if sp2 not in species[i]:
                                        weights[i] += 1
                                        species[i].append(sp2)
                                break
    if len(lists_contained) == 0:
        return lists_contained,[w*1./(len(in_species)+len(out_species)) for w in weights],species
    max_len = 0
    index = -1
    for i in range(len(lists_contained)):
        if len(lists_contained[i]) > max_len:
            index = i
            max_len = len(lists_contained[i])
    return [lists_contained[index]],[weights[index]*1./(len(in_species)+len(out_species))],[species[index]]            



# Count how many ingroup species is an interval preserved in.
# If each adjacency is preserved in at least one outgroup species, good.

def getCounts(block,in_species,out_species=[],to_print=False):
    count = {}
    association = {}
    out_association = {}
    for b in block:
        sp = b[:b.find('.')]
        if sp in in_species:
            if block[b] not in association.values():
                count[sp] = 1
                association[sp] = block[b]
                if to_print:
                    print b,association[sp]
            else:
                for k in association:
                    if association[k] == block[b]:
                        count[k] += 1
                        if to_print:
                            print b,association[k]
        if sp in out_species:
            out_association[sp] = block[b]
            if to_print:
                print b,out_association[sp]
    return count,association,out_association                        


