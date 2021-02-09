import itertools
import time

from halp.directed_hypergraph import DirectedHypergraph


#######################################################################################################################################################################

def find_cuts(H,taildistancelist):
    '''
    Finds the head and tail cuts based off the heuristic path's tail distance list.
    Returns the union of the head and tail cuts (represented by the edges which cross those cuts)
        the tail cuts in the same representation
        the head cuts in the same representation
    '''
    print('finding cuts')
    tailcuts = find_tail_cuts(H,taildistancelist) #normal cuts have the list of vertices in C and a distance value associated with them
    headcuts = find_head_cuts(H,taildistancelist)

    edgetailcuts = convert_vertex_cuts_to_edge_cuts(H,finalcuts)
    edgeheadcuts = convert_vertex_cuts_to_edge_cuts(H,headcutstuples)

    edgefinalcuts = edgetailcuts.union(edgeheadcuts)
    return edgefinalcuts,edgetailcuts,edgeheadcuts
            
def convert_vertex_cuts_to_edge_cuts(H,cutlist):
    '''
    Takes cuts in a cutlist that are defined by vertices in H and converts them into cuts that are based on edges in H
    '''
    edgecuts = set()
    for cut in range(len(cutlist)):
        E = find_crossing_edges(H,cut,cutlist)
        edgecuts.add(tuple(E))
    return edgecuts
        
def find_crossing_edges(H,cj,tailcuts):
    '''
    Finds the edges that cross the cut cj
    '''
    E = []
    #first make sure the entire tail is in cj
    Cj = tailcuts[cj][0]
    for e in H.hyperedge_id_iterator():
        tailin = True
        for v in H.get_hyperedge_tail(e):
            if v not in Cj:
                tailin = False
                break
        #Then make sure that at least one vertex in the head is not in cj
        if tailin == True:
            for v in H.get_hyperedge_head(e):
                if v not in Cj:
                    E.append(e)
                    break

    return E

def find_head_cuts(H,taildistancelist):
    '''
    Finds the head cuts which just union the vertices with smaller head distance for each head distance value.
    '''
    C = []
    V = set(['SUPERSOURCE'])
    #NOTE: This must be changed to accomodate weighted hypergraphs. Currently it uses 1 as the weight for each edge!
    distances = [taildistancelist[a] + 1 for a in taildistancelist]
    distances = sorted(set(distances))
    distances.append(1000000000000)
    #print('distances',distances)
    vertexdistancelist = {}
    for v in H.get_node_set():
        mindist = 100000
        for e in H.get_backward_star(v):
            if taildistancelist[e] + H.get_hyperedge_weight(e) < mindist:
                mindist = taildistancelist[e] + H.get_hyperedge_weight(e)
        vertexdistancelist[v] = mindist
            
    for d in distances:
        for v in H.get_node_set():
            if v in vertexdistancelist and vertexdistancelist[v] < d and v != 'SUPERTARGET':
                V.add(v)
        C.append(tuple(V))

    return C

def find_tail_cuts(H,taildistancelist):
    '''
    Finds the tail cuts which just union the tails of edges with the same or smaller tail distance for all tail distance values in the sorted tail distance list.
    '''
    C = []
    V = set(['SUPERSOURCE'])
    distances = [taildistancelist[a] for a in taildistancelist]
    distances = sorted(set(distances))
    distances.append(1000000000000)
    #print('distances',distances)
    for d in distances:
        for e in H.hyperedge_id_iterator():
            if e in taildistancelist and taildistancelist[e] < d:
                for v in H.get_hyperedge_tail(e):
                    V.add(v)
        C.append((tuple(V),d))

    return C


#######################################################################################################################################################################

