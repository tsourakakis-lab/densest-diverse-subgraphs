from collections import defaultdict
from FibonacciHeap import FibonacciHeap
import math
import numpy as np
from spectral_algorithms import *

def densest_subgraph_charikar_at_least_k( edges, k, original = True ):
     # 2 approximation - [Charikar, 2000]
    degrees = defaultdict(float)  # nodes and degrees
    adjlist = defaultdict(set)
    node_dict = dict()
    if original:
        for e in edges:
            u = int(e[0]); v = int(e[1]);
            w = 1
            degrees[u] += w
            degrees[v] += w
            adjlist[u].add(v)
            adjlist[v].add(u)
    else:
        for k,w in edges.items():
            u,v = k
            degrees[u] += w
            degrees[v] += w
            adjlist[u].add(v)
            adjlist[v].add(u)
    H = FibonacciHeap()
    for node in degrees:
        node_dict[node] = H.insert( degrees[node], node )
    total_weight = sum( degrees.values() )/2
    nodes = len(degrees)
    nodeSet = set([nd for nd in node_dict])
    maxSubgraph = nodeSet.copy()
    maxDensity = total_weight / nodes
    
    while H.n > 1 :
        if H.n <= k:
            break
        # peel node with lower (weighted) degree
        node = H.extract_min()
        nodeId, nodeDegree = node.value, node.key
        nodeSet.remove( nodeId )
        for neighbor in adjlist[nodeId]:
            if original: 
                w = 1
            else:
                w = edges[(min(nodeId,neighbor),max(nodeId,neighbor))]
            degrees[neighbor] -= w
            H.decrease_key(node_dict[neighbor], degrees[neighbor])
            adjlist[neighbor].remove( nodeId )
            total_weight -= w
        nodes -= 1
        if (total_weight / nodes) > maxDensity: 
            maxSubgraph = nodeSet.copy()
            maxDensity = total_weight / nodes
    return list(maxSubgraph), maxDensity
    

def Algorithm2( G, l, a ):
    C = max(l) + 1
    assert a >= 1/C
    edges = [e for e in G.edges()]
    k = math.ceil(1/a)
    S,dsalk_density = densest_subgraph_charikar_at_least_k( edges, k)
    return diversify2(set(S),l,G,a), dsalk_density


def diversify2(S,l,G,a):
    colors = defaultdict( set )
    color_sizes = defaultdict( int )
    # Init
    for c in np.unique(l):
        colors[c] = set()
        color_sizes[c] = 0
    for node in S:
        colors[l[node]].add(node)
        color_sizes[l[node]] += 1
    colors_outside = defaultdict( set )
    for node in G.nodes():
        if node not in S:
            colors_outside[l[node]].add(node)
    # Start by adding nodes outside subgraph
    try:
        while max(color_sizes.items(), key=lambda x:x[1])[1]/len(S) > a:
            # Drop colors that do not have any node outside S
            saturated_colors = [] # Keep track of all colors that we do not have any nodes outside
            for clr,clr_cnt in colors_outside.items():
                if len(clr_cnt) == 0:
                    saturated_colors.append( clr )
            for clr in saturated_colors:
                del color_sizes[clr]
            assert len(color_sizes) > 1 # We have saturated all but largest, throw error
            color = min(color_sizes, key=color_sizes.get)
            candidate_nodes = colors_outside[color]
            cn = [(v,len(set(G.neighbors(v)).intersection(S))) for v in candidate_nodes]
            new_node = max(cn,key=lambda item:item[1])[0]
            # Add new node to S
            S.add(new_node)
            colors_outside[color].remove(new_node)
            colors[color].add(new_node)
            color_sizes[color] += 1
    except:
        if True:
            # Rebuild counts
            color_sizes = defaultdict( int )
            for node in S:
                color_sizes[l[node]] += 1
            # If not possible to add more nodes, start pruning
            while max(color_sizes.items(), key=lambda x:x[1])[1]/len(S) > a:
                color = max(color_sizes, key=color_sizes.get)
                candidate_nodes = colors[color]
                cn = [(v,len(set(G.neighbors(v)).intersection(S))) for v in candidate_nodes]
                remove_node = min(cn,key=lambda item:item[1])[0]
                # Remove node to S
                S.remove(remove_node)
                colors[color].remove(remove_node)
                color_sizes[color] -= 1
        #except:
        #    print(colors)
        #    raise
    return S, colors

# Algorithms from Anagnostopoulos paper
def algo_2dfsg( G, l, t=0.01 ):
    # t is a tolerance parameter
    C = max(l)+1
    DS = flowless(G,5)[0]
    S, colors = diversify(DS,l,G,1/C+t)
    return S, colors

def power_iteration(A, num_simulations: int):
    # Ideally choose a random vector
    # To decrease the chance that our vector
    # Is orthogonal to the eigenvector
    b_k = np.random.rand(A.shape[1])

    for _ in range(num_simulations):
        # calculate the matrix-by-vector product Ab
        b_k1 = np.dot(A, b_k)

        # calculate the norm
        b_k1_norm = np.linalg.norm(b_k1)

        # re normalize the vector
        b_k = b_k1 / b_k1_norm

    return b_k

def single_sweep( A, G, l, fair=True, tol=0.01 ):
    n = A.shape[0]
    n_colors = max(l)+1
    #v = power_iteration( M, 1000)
    map_node_color = returnIndexColorMap( l )
    if fair:
        if max(l) > 1: # More than two colors
            v = compute_main_eigenvector_fairified_matrix__scipy__MORE_THAN_TWO_COLORS(A, map_node_color)[0]
        else:
            v = compute_main_eigenvector_fairified_matrix__scipy(A, map_node_color)[0]
    else:
        v = compute_main_eigenvector__scipy( A )[0]
    v = np.array(v[0].tolist())
    S, tempS, DS, newDS = set(), set(), 0, 0
    # Non-increasing
    sorted_v = np.argsort(-v)
    color_count = defaultdict( int )
    for s in range(n):
        new_node = sorted_v[s]
        color_count[l[new_node]] += 1
        tempS.add(new_node)
        newDS = (newDS + len(set(G.neighbors(new_node)).intersection(tempS))) / (s+1)
        if ( newDS > DS ) and ( ( max(color_count.values()) / (s+1) ) <= 1/n_colors + tol ):
            DS = newDS
            S = tempS.copy()
    # Non-decreasing
    tempS, newDS = set(), 0
    sorted_v = np.argsort(v)
    color_count = defaultdict( int )
    for s in range(n):
        new_node = sorted_v[s]
        color_count[l[new_node]] += 1
        tempS.add(new_node)
        newDS = (newDS + len(set(G.neighbors(new_node)).intersection(tempS))) / (s+1)
        if ( newDS > DS ) and ( ( max(color_count.values()) / (s+1) ) <= 1/n_colors + tol ):
            DS = newDS
    # Non-increasing absolute value
    tempS, newDS = set(), 0
    sorted_v = np.argsort(-np.absolute(v))
    color_count = defaultdict( int )
    for s in range(n):
        new_node = sorted_v[s]
        color_count[l[new_node]] += 1
        tempS.add(new_node)
        newDS = (newDS + len(set(G.neighbors(new_node)).intersection(tempS))) / (s+1)
        if ( newDS > DS ) and ( ( max(color_count.values()) / (s+1) ) <= 1/n_colors + tol ):
            DS = newDS
    # Non-decreasing absolute value
    tempS, newDS = set(), 0
    sorted_v = np.argsort(np.absolute(v))
    color_count = defaultdict( int )
    for s in range(n):
        new_node = sorted_v[s]
        color_count[l[new_node]] += 1
        tempS.add(new_node)
        newDS = (newDS + len(set(G.neighbors(new_node)).intersection(tempS))) / (s+1)
        if ( newDS > DS ) and ( ( max(color_count.values()) / (s+1) ) <= 1/n_colors + tol ):
            DS = newDS
    return S,DS

def pairedSweep( A, G, l, fair=True ):
    n = A.shape[0]
    n_colors = max(l) + 1
    #v = power_iteration( M, 1000)
    map_node_color = returnIndexColorMap( l )
    if fair:
        if max(l) > 1: # More than two colors
            v = compute_main_eigenvector_fairified_matrix__scipy__MORE_THAN_TWO_COLORS(A, map_node_color)[0]
        else:
            v = compute_main_eigenvector_fairified_matrix__scipy(A, map_node_color)[0]
    else:
        v = compute_main_eigenvector__scipy( A )[0]
    v = np.array(v[0].tolist())
    colors = defaultdict( set )
    color_sizes = defaultdict( int )
    for node in G.nodes():
        colors[l[node]].add(node)
        color_sizes[l[node]] += 1
    min_c = min(color_sizes.values())
    # Sort every vector by each color
    S, tempS, DS, newDS = set(), set(), 0, 0
    sorted_v = np.argsort(-v)
    color_split = defaultdict(list)
    for s in range(n):
        s_color = l[sorted_v[s]]
        color_split[s_color].append(sorted_v[s])
        
    for s in range(min_c):
        t = 0
        for c in color_split:
            new_node = color_split[c][s]
            t += len(set(G.neighbors(new_node)).intersection(tempS))
            tempS.add( new_node )
        newDS += t
        if newDS/(n_colors*(s+1)) > DS:
            DS = newDS/(n_colors*(s+1))
            S = tempS.copy()
    
    # Non-decreasing
    tempS, newDS = set(), 0
    sorted_v = np.argsort(v)
    color_split = defaultdict(list)
    for s in range(n):
        s_color = l[sorted_v[s]]
        color_split[s_color].append(sorted_v[s])
    for s in range(min_c):
        t = 0
        for c in color_split:
            new_node = color_split[c][s]
            t += len(set(G.neighbors(new_node)).intersection(tempS))
            tempS.add( new_node )
        newDS += t
        if newDS/(n_colors*(s+1)) > DS:
            DS = newDS/(n_colors*(s+1))
            S = tempS.copy()
    # Non-increasing absolute value
    tempS, newDS = set(), 0
    sorted_v = np.argsort(-np.absolute(v))
    color_split = defaultdict(list)
    for s in range(n):
        s_color = l[sorted_v[s]]
        color_split[s_color].append(sorted_v[s])
    for s in range(min_c):
        t = 0
        for c in color_split:
            new_node = color_split[c][s]
            t += len(set(G.neighbors(new_node)).intersection(tempS))
            tempS.add( new_node )
        newDS += t
        if newDS/(n_colors*(s+1)) > DS:
            DS = newDS/(n_colors*(s+1))
            S = tempS.copy()
    # Non-decreasing absolute value
    tempS, newDS = set(), 0
    sorted_v = np.argsort(np.absolute(v))
    color_split = defaultdict(list)
    for s in range(n):
        s_color = l[sorted_v[s]]
        color_split[s_color].append(sorted_v[s])
    for s in range(min_c):
        t = 0
        for c in color_split:
            new_node = color_split[c][s]
            t += len(set(G.neighbors(new_node)).intersection(tempS))
            tempS.add( new_node )
        newDS += t
        if newDS/(n_colors*(s+1)) > DS:
            DS = newDS/(n_colors*(s+1))
            S = tempS.copy()
            
    return S, DS

def returnIndexColorMap( labels ):
    map__index__color = {}
    for i in range(len(labels)):
        map__index__color[i] = labels[i]
    return map__index__color
            
        
    

if __name__ == '__main__':
    pass