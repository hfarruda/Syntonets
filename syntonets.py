#Codes created by Henrique Ferraz de Arruda and Luciano da Fontoura Costa.
#F. Costa, Luciano, and F. de Arruda, Henrique. "Syntonets: Toward A Harmony-Inspired General Model of Complex Networks." arXiv preprint arXiv:1910.11047v2 (2019).
import numpy as np 
import igraph as ig
import os 
import sys
import random

def _helmholtz_test(x,y,f1,f_delta_min=2,f_delta_max=38):
    fX = np.array(list(x.keys()),dtype = np.float)
    MX = np.array(list(x.values()),dtype = np.float)
    fY = np.array(list(y.keys()),dtype = np.float)
    MY = np.array(list(y.values()),dtype = np.float)
    #Comparing all pairs
    d_i = []
    c_i = []
    for i,f in enumerate(fX):
        diff = np.abs(fY-f)
        consonance_pos = np.argwhere((diff < f_delta_min/2.) == True)
        dissonance_pos = np.argwhere(((diff >= f_delta_min/2.) * (diff < f_delta_max/2.)) == 1)
        
        #Computing conssonance
        mx = MX[i]
        if consonance_pos.shape[0] != 0:
            consonance_pos = consonance_pos.T[0]
            aux = 0
            for p in consonance_pos:
                aux += mx * MY[p]
            delta = 1
            c_i.append(aux * delta)
        else:
            c_i.append(0.)
        
        #Computing dossonance
        if dissonance_pos.shape[0] != 0:
            dissonance_pos = dissonance_pos.T[0]
            aux = 0
            for p in dissonance_pos:
                aux += mx * MY[p]
            delta = 1
            d_i.append(aux * delta)
        else:
            d_i.append(0.)

    C = np.sum(c_i)
    D = np.sum(d_i)
    return {'D':D, 'C':C}


def _create_freq2M_dict(f1, beta, alpha, min_frequency=20, max_frequency=20000):
    f_1h_i = []
    h = []
    i = 1
    aux = f1 * np.power(float(i), beta)
    while (aux<=max_frequency and aux>=min_frequency):
        f_1h_i.append(aux)
        h.append(i)
        i += 1
        aux = f1 * np.power(float(i), beta)
    M = np.exp(-alpha * np.array(h))
    f2M = {f_1h_i[i]:M[i] for i in range(len(M))}
    return f2M

def _adaptive_threshold(g, giant, number_of_edges):
    weights = g.es["weight"]
    sortedWeights = np.sort(weights)[::-1]
    if number_of_edges > 0:
        if number_of_edges >= len(sortedWeights):
            number_of_edges = len(sortedWeights) - 1
        threshold = sortedWeights[number_of_edges]
        delete_edges = [i for i,w in enumerate(weights) if w <= threshold]
        g.delete_edges(delete_edges)
    
    if giant:
        g = g.components(mode=ig.WEAK).giant()

    return g

def create_network(interval2note, interval2frequency, beta = 1, alpha = 0.2, number_of_edges = -1, giant = True, with_colors = True):
    """
    Create the consonance-based network.
    For more information see: https://arxiv.org/pdf/1910.11047.pdf

    Parameters
    -------
    interval2note : list representing the notes respectively to the list positions.
    Input : list of strings

    interval2frequency : list representing the frequencies respectively to the list
    positions.
    Input : list of float numbers

    beta : parameter of the consonance network that induces a partial shift.
    The default value is beta = 1, which recovers the perfectly 
    harmonic partials.
    Input: float

    alpha : parameter of the consonance network that controls the exponential
    decay. The default value is alpha = 0.2.
    Input : float
    
    number_of_edges : a parameter that controls the number of edges in the 
    resultant network. The default value is beta = -1 that keeps all the 
    original edges.
    Input : integer
    
    giant : if True, only the weakly connected component is kept. The default 
    value is giant = True. In the case of True, the number of edges can vary.
    Input : bool
    
    with_colors : if True, collors are assigned as node propoerties (e.g., 
    g.vs['color']). These colors are selected for each of the 
    note.
    Input : bool

    Returns
    -------
    (Igraph graph with the consonance network,
     Igraph graph with the dissonance network)
    """
    intervals = [i for i in range(len(interval2note))]
    f1 = interval2frequency[0]
    note2frequencies = {interval2note[i]:_create_freq2M_dict(interval2frequency[i], beta, alpha) for i in intervals}
    
    D = dict()
    C = dict()
    for pos,i in enumerate(note2frequencies.keys()):
        for j in list(note2frequencies.keys())[pos+1::]:
            x = note2frequencies[i]
            y = note2frequencies[j]
            D_and_C = _helmholtz_test(x,y,f1)
            D[(i,j)] = D_and_C["D"]
            C[(i,j)] = D_and_C["C"]
    
    #With Igraph
    g_C = ig.Graph()
    g_C.add_vertices(interval2note)
    
    g_D = ig.Graph()
    g_D.add_vertices(interval2note)

    if with_colors:
        colors = ["#e93f33", "#f3a433", "#a52a2a", "#fdf731", "#fefde0", "#7ff21a", "#90ee90", "#76f9fb", "#add8e6", "#244afc", "#e981ee", "#f4bfca"]
        g_C.vs['color'] = [colors[i%12] for i in range(g_C.vcount())]
        g_D.vs['color'] = [colors[i%12] for i in range(g_D.vcount())]
    
    edges = list(C.keys())
    weights = list([C[edge] for edge in edges])
    g_C.add_edges(edges) 
    g_C.es["weight"] = weights

    edges = list(D.keys())
    weights = list([D[edge] for edge in edges])
    g_D.add_edges(edges) 
    g_D.es["weight"] = weights

    g_C = _adaptive_threshold(g_C, giant, number_of_edges)
    g_D = _adaptive_threshold(g_D, giant, number_of_edges)
    
    return g_C, g_D


def visualize(g, fig_name, vertex_size = 40, edge_width_scale = -1, niter = 1000, with_colors = True):
    """
    Visuliaze the consonance-based network by employing Igraph library.

    Parameters
    -------
    g : graph that represents the consonance-based network.
    Input: igraph graph

    fig_name : set the name of the output file.
    Input: string

    vertex_size : defines the node size of the visualization. The default 
    value is vertex_size = 40.
    Input:

    edge_width_scale : scale the edge width according to their weights.
    Input: integer

    niter : defines the number of iterations of the Fruchterman Reingold 
    function. The default value is niter = 1000.
    Input: integer

    with_colors : if with_colors = True, the visualization uses
    the node attribute "color" defened in the graph (g.vs["color"]).
    Input: bool

    Output
    -------
    Save a 'png' figure.
    """
    layout = g.layout_fruchterman_reingold(weights=g.es["weight"], niter = niter)
    visual_style = {}
    visual_style["vertex_size"] = vertex_size
    if with_colors:
        visual_style["vertex_color"] = g.vs["color"]
    visual_style["vertex_label"] = g.vs["name"]
    if edge_width_scale > 0:
        visual_style["edge_width"] = [edge_width_scale * w for w in g.es["weight"]]
    visual_style["layout"] = layout
    visual_style["bbox"] = (800, 800)
    visual_style["margin"] = 30

    ig.plot(g, target = fig_name + ".png", **visual_style)