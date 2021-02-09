from xml.dom import minidom
import cplex

from halp.directed_hypergraph import DirectedHypergraph
from halp.algorithms.directed_paths import b_visit

def make_shortest_cyclic_hyperpath_ilp(H,source,target,outfile):
    out = open(outfile,'w')

    ## write objective
    out = writeObjective(H,out)

    ## write constraints
    out.write('Subject To\n')

    #THIS IS WHAT MAKES IT THE ILP

    #out = writeBinaryBounds(H,out)

    out.write('End\n')
    out.close()


def writeObjective(H,out,minimize=True):
    if minimize:
        out.write('Minimize\n')
    else:
        out.write('Maximize\n')

    for hedge in H.hyperedge_id_iterator():
        out.write(' + %d %s' % (H.get_hyperedge_attribute(hedge,'weight'),a(hedge)))
    out.write('\n')
    return out

def writeBinaryBounds(H,out,write_q=False,write_f=None):
    '''
    Specify all the alpha variables as binary.
    '''
    out.write('Binary\n')
    for hnode in H.node_iterator():
        out.write(' %s\n' % (a(hnode)))
    for hedge in H.hyperedge_id_iterator(): # for all e \in E
        out.write(' %s\n' % (a(hedge)))
        if write_q:
            out.write(' %s\n' % (q(hedge)))
    if write_f:
        ## write_f is a GRAPH
        edges = [(write_f.get_hyperedge_tail(e)[0],write_f.get_hyperedge_head(e)[0]) for e in write_f.get_hyperedge_id_set()]
        for e in edges:
            out.write(' %s\n' % (f(e[0],e[1])))
    return out


def t_in_S(S,target):
    #here is some code about checking if the sink is in S
    #returns a boolean
    for node in S:
        if node == target:
            return True
    return False

def get_new_constraint(S,H):
    #this takes the S,T cut and the set of edges and finds all edges that cross the cut
    #this one is a little harder to think through, but this is just the exponential version.  It should be fine for our hypergraphs, but does not work for any hypergraph.
    #returns a list of edges that cross the s,t-cut
    crossedges = []
    addedge = False
    nodesintail = True
    for hedge in H.hyperedge_id_iterator():
        addedge = False
        nodesintail = True
        for tailnode in H.get_hyperedge_tail(hedge):
            if tailnode not in S:
                nodesintail = False
        #if it got here, that means that all the nodes in the tail are in S
        if nodesintail == True:
            for headnode in H.get_hyperedge_head(hedge):
                if headnode not in S:
                    #this means that this edge crosses the cut
                    addedge = True
            if addedge == True:
                crossedges.append(hedge)
    return crossedges


def reachability_from_edges(H,ones,source):
    #takes the solution from the previous iteration of the ILP and sees what nodes are reachable from those edges.  Make a subhypergraph, and then call b_visit on it.
    #returns a list of nodes, S
    #print( '_____________rfe_______')
    #print( ones)
    nodeset = set()
    H2 = DirectedHypergraph()
    H2.add_nodes(H.get_node_set())
    for edge in ones:
        H2.add_hyperedge(H.get_hyperedge_tail(una(edge)),H.get_hyperedge_head(una(edge)))
    
    nodeset = b_visit(H2,source)
    #print ('nodeset size:',len(nodeset[0]), nodeset[0])
    return nodeset


def get_addable_cuts(cuts):
    '''

    '''
    constraints = []
    for cut in cuts:
        constraint = make_constraint_from_edges(cut)
        constraints.append(constraint)
    return constraints

def make_constraint_from_edges(crossedges):
    '''

    '''

    #addconstraint to the ILP
    crossedges = set(crossedges)
    ones = [a(c) for c in crossedges]
    val = [1.0]*len(ones)
    eq = cplex.SparsePair(ind=ones,val=val)
    return eq


def ILP_from_LP(lp,H,nodeset):
    '''

    '''
    #print('-'*30,'ILP')
    for e in H.hyperedge_id_iterator():
        lp.variables.set_types(a(e),lp.variables.type.binary)
    lp.solve()
    #print('\n'*30)
    print('ILP lower bound is {}'.format(lp.solution.get_objective_value()))
    #for e in H.hyperedge_id_iterator():
        #if lp.solution.get_values(a(e)) == 1:
            #print(e,H.get_hyperedge_tail(e),H.get_hyperedge_head(e))

    #print('\n'*30)
    return lp

def build_initial_LP(lpfile,cuts):
    '''

    '''
    lp = cplex.Cplex()
    lp.read(lpfile)
    I = get_addable_cuts(cuts)
    #print(I)
    for j in I:
        lp.linear_constraints.add(lin_expr = [j], senses = ['G'], rhs = [1], names = ['starter{0}'.format(j)])
    #print ('______________numconstraints____________')
    #print(len(I))
    return lp


def run_ILP_cuttingplanes(H,nodeset,lpfile,outprefix,numsols,target,source,cuts,normalcuts,headcuts,targetname='SUPERTARGET'):
    '''

    '''
    #print(cuts)
    edgedict = {}
    for e in H.hyperedge_id_iterator():
        edgedict[e] = 0

    lp = build_initial_LP(lpfile,headcuts.union(normalcuts))
    lp.set_results_stream(logfile)
    lp.set_error_stream(logfile)
    S = [source]
    ilp = ILP_from_LP(lp,H,nodeset)

    ilp.set_warning_stream(logfile)
    numiterations = solveILPcuttingplanes(H,ilp,target,source,outprefix,nodeset,numitrs = 25000,targetname=targetname,ilpchar='b')

    
def solveILPcuttingplanes(H,ilp,target,source,outprefix,nodeset,numitrs=100000,verbose=False,targetname='SUPERTARGET',ilpchar='2'):
    '''

    '''
    numsolsfound = 1
    numoptobjective = 0
    maxobj = None
    allvars = []
    S = [source]
    numiterations = 0
    ofile = open('{}ilpdata{}.txt'.format(targetname[-6:],ilpchar),'w')
    while t_in_S(S,target) == False and numiterations < numitrs:
        crossedges = get_new_constraint(S,H)
        ones = [a(c) for c in crossedges]
        val = [1.0]*len(ones)
        eq = cplex.SparsePair(ind=ones,val=val)
        ilp.linear_constraints.add(lin_expr = [eq], senses = ['G'], rhs = [1], names = ['iteration%d' % (numiterations)])
        ilp.solve()
        
        if ilp.solution.pool.get_num()>0:
            objective = ilp.solution.pool.get_objective_value(0)
            ofile.write('iteration {} solval {}\n'.format(numiterations,objective))
            ilp.solution.pool.write('%s-%d.sol' % (outprefix,numsolsfound),0)
            variables = getILPSolution(H,nodeset,outprefix,numsolsfound,objective,verbose)
            ones = [var for var in variables if variables[var] == 1 and var[0:3]=='a_e' and var.split('_')[1] not in nodeset]
            S,temp1,temp2,temp3 = reachability_from_edges(H,ones,source)
            S.add(source)
        else:
            print ('Infeasible Solution. quitting.')
            break
        numsolsfound+=1
        numiterations+=1

    allvars = variables
    #print ('-'*20 + 'Cplex Output End' + '-'*20 + '\n')
    #print ('%d solutions found' % (numsolsfound-1))
    return numiterations


def getILPSolution(H,nodeset,outprefix,num,objective,verbose):
    #print ('\nGetting ILP Solution for Solution # %d in Pool' % (num))

    # parse xml
    xml = minidom.parse('%s-%d.sol' % (outprefix,num))
    cplexsol = xml.getElementsByTagName('CPLEXSolution')[0]
   
    # get variables
    elements = cplexsol.getElementsByTagName('variables')[0].getElementsByTagName('variable')
    variables = {}
    for v in elements:
        variables[str(v.getAttribute('name'))] = float(v.getAttribute('value'))

    out = open('%s-%d.variables' % (outprefix,num),'w')
    out.write('# Objective = %s\n' % (objective))
    out.write('#name\tval\trounded_val\n')

    if verbose:
        print ('VARIABLES:')
    
    numrounded = 0
    for v in variables:
        rounded = int(variables[v]+0.5)
        if verbose and variables[v] != 0:
            print(v,variables[v],rounded)
        out.write('%s\t%f\t%d\n' % (v,variables[v],rounded))
    out.close()
    
    #print (' wrote variables to file %s' % ('%s-%d.variables' % (outprefix,num)))

    return variables
    

def a(name):
    return 'a_%s' % (name)

def q(hedge):
    return 'q_%s' % (hedge)

def f(node1,node2):
    return 'f_%s_%s' % (node1,node2)

def una(name):
    return name[2:].strip()
