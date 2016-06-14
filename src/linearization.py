import sys
import os
import math
import networkx as nx


# Create auxiliary graph based on mults
def buildAuxGraph(G,multiplicity,f):
    auxG=nx.Graph()
    # Step 1: Adding 2 vertices for each copy of a marker
    for v in auxG:
        for k in range(0,f*multiplicity[v]):
            auxG.add_node(str(v)+'_'+str(k))
    # Step 2: Adding 2 vertices for each adjacency.
    # Also adding corresponding edges and giving them weights as given in the matrix
    edge_found = {}
    for i in range(len(G.edges())):
        e = G.edges()[i]
        #print "Edges",e
        edge_found[e] = 0 
        #print G[e[0]][e[1]]['weight']
        if e[0] == e[1]:
            for x in e:
                auxG.add_node('e'+str(i)+'_'+str(x)+'_tandem0')
                auxG.add_node('e'+str(i)+'_'+str(x)+'_tandem1')
                for k in range(0,f*multiplicity[x]):
                    auxG.add_edge(str(x)+'_'+str(k),'e'+str(i)+'_'+str(x)+'_tandem0')
                    auxG.add_edge(str(x)+'_'+str(k),'e'+str(i)+'_'+str(x)+'_tandem1')
                    auxG[str(x)+'_'+str(k)]['e'+str(i)+'_'+str(x)+'_tandem0']['correspondence'] = e 
                    auxG[str(x)+'_'+str(k)]['e'+str(i)+'_'+str(x)+'_tandem1']['correspondence'] = e 
                    auxG[str(x)+'_'+str(k)]['e'+str(i)+'_'+str(x)+'_tandem0']['weight'] = G[e[0]][e[1]]['weight']
                    auxG[str(x)+'_'+str(k)]['e'+str(i)+'_'+str(x)+'_tandem1']['weight'] = G[e[0]][e[1]]['weight']
                auxG.add_edge('e'+str(i)+'_'+str(e[0])+'_tandem0','e'+str(i)+'_'+str(e[1])+'_tandem1')
                auxG['e'+str(i)+'_'+str(e[0])+'_tandem0']['e'+str(i)+'_'+str(e[1])+'_tandem1']['correspondence'] = [] 
                auxG['e'+str(i)+'_'+str(e[0])+'_tandem0']['e'+str(i)+'_'+str(e[1])+'_tandem1']['weight'] = G[e[0]][e[1]]['weight']
            continue
        for x in e:
            auxG.add_node('e'+str(i)+'_'+str(x))
            for k in range(0,f*multiplicity[x]):
                auxG.add_edge(str(x)+'_'+str(k),'e'+str(i)+'_'+str(x))
                auxG[str(x)+'_'+str(k)]['e'+str(i)+'_'+str(x)]['correspondence'] = e 
                auxG[str(x)+'_'+str(k)]['e'+str(i)+'_'+str(x)]['weight'] = G[e[0]][e[1]]['weight']
        auxG.add_edge('e'+str(i)+'_'+str(e[0]),'e'+str(i)+'_'+str(e[1]))
        auxG['e'+str(i)+'_'+str(e[0])]['e'+str(i)+'_'+str(e[1])]['correspondence'] = [] 
        auxG['e'+str(i)+'_'+str(e[0])]['e'+str(i)+'_'+str(e[1])]['weight'] = G[e[0]][e[1]]['weight']
    return auxG,edge_found


def linearization(G,multiplicity,f): # multiplicity is a dictionary, with markers and their reversals as keys.
    # Creating auxiliary adjacency graph
    auxG,edge_found = buildAuxGraph(G,multiplicity,f)
    #for e in auxG.edges():        
    #    print e[0],e[1],auxG[e[0]][e[1]]
    # Maximum weight matching algorithm called on the auxiliary graph.
    # This gives us A maximum matching, and, correspondingly, a maximum weight 2m-matching in the matrix
    # This DOES NOT maximize the number of adjacencies chosen.
#    print "Here",auxG.edges()
    matching=nx.max_weight_matching(auxG)
#    print "Done",matching
    present=0

    # Using the matching to find a component mCi1P respecting spanning subgraph of maximum weight.
    # The corresponding matrices- for the adjacencies retained, and those discarded- are the outputs.
    to_keep = []
    #traversed = []
    #for v in matching:
    #    #if v not in traversed and matching[v] not in traversed:
    #    #    print v,matching[v]
    #    #    traversed.append(v)
    #    #    traversed.append(matching[v])
    for vertex in matching.keys():
        #print vertex,matching[vertex],auxG[vertex][matching[vertex]]['correspondence']
        if not (vertex[0]=='e'):
            # Condition for retaining an adjacency: both ends of the corresponding adjacency must be matched, but not to each other.
            edge_id = ''
            marker_id = ''
            edge_name = ''
            edge_found[auxG[vertex][matching[vertex]]['correspondence']] += 1
#            print vertex,matching[vertex],auxG[vertex][matching[vertex]]['correspondence'],edge_found[auxG[vertex][matching[vertex]]['correspondence']]
            if edge_found[auxG[vertex][matching[vertex]]['correspondence']] >= 2 and not(auxG[vertex][matching[vertex]]['correspondence'] in to_keep):
                to_keep.append(auxG[vertex][matching[vertex]]['correspondence'])
#            if vertex[0] == 'e':
#                edge_id = vertex[:2]
#                edge_name = vertex
#                marker_id = matching[vertex][:matching[vertex].find('_')]
#            else:
#                edge_id = matching[vertex][:2]
#                edge_name = vertex
#                marker_id = vertex[:vertex.find('_')]
#            e = auxG[vertex][matching[vertex]]['correspondence']
#            for v in matching.keys():
#                if v[:2] == edge_id and v != edge_name:
#                    if auxG[v][matching[v]]['correspondence'] != None:
#                        if not(auxG[vertex][matching[vertex]]['correspondence'] in to_keep):
#
#            for e in edge_itemizing:
#                if u == edge_itemizing
##                print 
#                to_keep.append(auxG[vertex][matching[vertex]]['correspondence'])
#            first_end=None
#            other_end=None
#            if vertex[0:5]=='e':
#                k=5
#                while vertex[k]!='_':
#                    k+=1
#                row=int(vertex[5:k])
#                first_end=int(vertex[k+1:])
#                            for col in m.get_row_info(row)._set:
#                                    if col!=first_end:
#                        other_end=col
#            elif matching[vertex][0:5]=='e':
#                k=5
#                while matching[vertex][k]!='_':
#                    k+=1
#                row=int(matching[vertex][5:k])
#
#                first_end=int(matching[vertex][k+1:])
#
#                for col in m.get_row_info(row)._set:
#                    if col!=first_end:
#                        other_end=col
#            # Check if both ends are matched
#            if  'e'+str(row)+'_'+str(other_end) in matching.keys():
#                if m.get_row_info(row) not in new_matrix:
#                    new_matrix.add_row_info(m.get_row_info(row))
#            else:
#            # Case when one edge-vertex is matched to a copy of a column, but the other is unmatched.
#                if m.get_row_info(row) not in discarded_rows:
#                    discarded_rows.add_row_info(m.get_row_info(row))
#        elif (vertex[0:5]=='e' and  matching[vertex][0:5]=='e'):
#            k=5
#        # The two edge-vertices of an adjacency are matched to each other.
#            while vertex[k]!='_':
#                k+=1
#            if m.get_row_info(int(vertex[5:k])) not in discarded_rows:
#                discarded_rows.add_row_info(m.get_row_info(int(vertex[5:k])))
#     
    return to_keep #[new_matrix,discarded_rows]
    
