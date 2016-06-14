#!/usr/bin/env


# ANGES 1.0, reconstruction ANcestral GEnomeS maps
# May 2012.
# Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

import sys
import os

import DataStructures.bm as bm

from numpy import *		# http://numpy.scipy.org/
from numpy.linalg import *
from decimal import *

#######################################################
#     ACS_weight_linear_interpolation
#
#     Given a species tree where a node is marked as the ancestor, 
#     computes the weight of a set of ACS (matrix) using linear interpolation (Fitch)
#
#######################################################

# recusive method for make_matrix
# buff - str, i - int, par - int, st - int, end - int,
# M - array of double, sp list of str
# return - int, int, int, int
def make_matrix_rec(buff, i, par, st, end, M, sp):
	curr = -1	# index of the current node
	A = - 1		# index in the matrix of the ancestor	
	if buff[i] == '(':		# internal node
		isLeaf = False		# True if the current node is a leaf	
		end -= 1		
		curr = end
		A, i, st, end = make_matrix_rec(buff, i + 1, curr, st, end, M, sp)
	elif buff[i] != ';':		# leaf
		isLeaf = True
		curr = st
		st += 1
		j = buff.index(':', i)		# index of next node
		sp.append(buff[i:j])
		i = j + 1
	# ancestor node	
	if i < len(buff) and buff[i] == '@':
		A = curr				
		i += 1	
	if not isLeaf:
		try:
			i = buff.index(':', i) + 1
		except ValueError:
			i = i+1
	if i < len(buff) and buff[i] != ';' and par >= 0:
		if isLeaf:
			M[curr, curr] = 1
		# find next node
		try:
			n = buff.index(',', i)		# index of next , node
		except ValueError:
			n =len(buff)
		next = min(n, buff.index(')', i))		# index of next node			
		# record length in matrix
		if par >= 0:
			l = float64(buff[i:next])		# length of edge			
			# if l is zero ignore the row
			if l == 0:
				# move M[curr]'s values to M[par]
				for ind in xrange(curr):
					M[par, ind] += M[curr, ind]
					M[ind, par] += M[ind, curr]				
				M[par, par] += M[curr, curr]
				# make curr apply the identity
				for ind in xrange(curr):
					M[curr, ind] = 0
					M[ind, curr] = 0
				M[curr, curr] = 1
			else:
				M[par, curr] = -1 / l
				M[par, par] += 1 / l				
				if not isLeaf:
					M[curr, curr] += 1 / l
					M[curr, par] = -1 / l
				
		if n == next:
			a, i, st, end = make_matrix_rec(buff, next + 1, par, st, end, M, sp)			
			if a >= 0:
				A = a
		else:
			i = next + 1	
	return A, i, st, end
#enddef

# makes the weight matrix of the tree in treeFile
# treeFile - str, offset - int
# return - array of double, list of str
def make_matrix(treeFile, offset):
	sp = []		# list of species in the tree
	f = file(treeFile, 'r')		# tree file=	
	buff = f.read()		# tree file contents
	f.close()
	n = buff.count(':') + 1		# number of nodes in the tree
	M = array([[float64(0.0) for j in xrange(n)] for i in xrange(n)])		# weight matrix
	A, i, st, end = make_matrix_rec(buff, 0, -1, 0, n,  M, sp)		# auxilary variables
	M = inv(M)
	return M[A - offset, :], sp
#enddef

# weights the martix m, according to the pylogenic tree in the file, treeFile
# set double to True if the matrix has doubled markers so that signle marker 
# rows are given 'infinite' weight
# m - bm.BinaryMatrix, treeFile - str, double - bool, offset - int
# return - bm.BinaryMatrix
def weight_matrix(m, treeFile, double, offset):
	# make inverse weight matrix from pylogenic tree
	M, sp = make_matrix(treeFile, offset)		# weight matrix, species in tree
	# calculate weight of each row
	for r in xrange(m._height):
		row = m.get_row_info(r)		# current row
		row._weight = Decimal('0')
		# add weight for each species in the row
		if row._isT or (not double or len(row._set) != 2 or len(set([(x + 1) // 2 for x in row._set])) > 1):
			for s in row._sp:
				try:		
					row._weight = row._weight + Decimal(str(M[sp.index(s)]))
				except ValueError:
					pass
		else:
			row._weight = Decimal('10000')
	return m
#enddef

def weight_intervals_local(intervals_list,species_dict,tree_file):
    M,sp = make_matrix(tree_file,0)
    weights = []
    for i in range(len(intervals_list)):
        weights.append(Decimal('0'))
        for s in species_dict[intervals_list[i]]:
            weights[i] += Decimal(str(M[sp.index(s)]))
    floated_weights = [float(str(x)) for x in weights]            
    return floated_weights
# intervals_list is a list of tuples, with intervals as first item 
# and list/tuple of species as the second

def weight_intervals(intervals_list,tree_file):
    M,sp = make_matrix(tree_file,0)
    weights = []
    for i in range(len(intervals_list)):
        weights.append(Decimal('0'))
        for s in intervals_list[i][2]:
            weights[i] += Decimal(str(M[sp.index(s)]))
    floated_weights = [float(str(x)) for x in weights]            
    return floated_weights
