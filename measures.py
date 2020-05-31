import networkx as nx
import statistics as stat
from sklearn.metrics.cluster import normalized_mutual_info_score


def density(graph, node2com, com2nodes, weight='weight'):
    """
    Formula:
        return_value = 1.0 / m * SIGMA { ||E(C_i)|| }, whereas:
        m : number of links in the graph
        E(C_i) : number of links that are inside community C_i

    Arguments:
        graph {[nx.Graph]} -- [network to be processed]
        node2com {[dict]} -- [keys are nodes, values are corresponding community index]
        com2nodes {[dict]} -- [keys are community indexes, values are nodes inside each community]
    
    Keyword Arguments:
        weight {str} -- [either 'weight' or 'prob'] (default: {'weight'})
    
    Returns:
        [float] -- [the density of network based on given partition]
    """
    sigma = 0
    m = graph.size(weight=weight)
    for edge in graph.edges(data=True):
        if node2com[edge[0]] == node2com[edge[1]]:
            sigma += edge.get(weight, 1.0)

    return float(sigma / m)


def internal_density(graph, node2com, com2nodes, ci, weight='weight'):
    """
    Formula:
        return_value = (1.0 / ||C_i||) * SIGMA{ mi / C(ni, 2)}, whereas:
        ||C_i|| = ni : number of nodes in community ci
        mi : number of links in community ci
        C(ni, 2) = ni * (ni - 1) / 2

    Arguments:
        graph {[nx.Graph]} -- [network to be processed]
        node2com {[dict]} -- [keys are nodes, values are corresponding community index]
        com2nodes {[dict]} -- [keys are community indexed, values are nodes inside each community]
        ci {[int]} -- [index of the intended community]
    
    Keyword Arguments:
        weight {str} -- ['weight' or 'prob'] (default: {'weight'})
    
    Returns:
        [float] -- [internal density of community ci of the graph based on given partition]
    """
    mi = dict()
    for edge in graph.edges(data=True):
        com = node2com[edge[0]]
        if node2com[edge[1]] == com:
            mi[com] = mi.get(com, 0) + edge.get(weight, 1.0)
    
    id_measure = 0
    for com, nodes in com2nodes.items():
        ni = len(nodes)
        id_measure += mi[com] / (ni * (ni - 1) / 2.0)
    
    return float(id_measure / len(com2nodes[ci]))


def edges_inside(graph, node2com, com2nodes, ci, weight='weight'):
    """
    Formula:
        return_value = (1.0 / ||C_i||) * SIGMA{ mi }, whereas:
        ||C_i|| : number of nodes in community ci
        mi : number of links in community ci

    Arguments:
        graph {[nx.Graph]} -- [network to be processed]
        node2com {[dict]} -- [keys are nodes, values are corresponding community index]
        com2nodes {[dict]} -- [keys are community indexed, values are nodes inside each community]
        ci {[int]} -- [index of the intended community]
    
    Keyword Arguments:
        weight {str} -- ['weight' or 'prob'] (default: {'weight'})
    
    Returns:
        [float] -- [edges inside of community ci of the graph based on given partition]
    """
    mi = dict()
    for edge in graph.edges(data=True):
        com = node2com[edge[0]]
        if node2com[edge[1]] == com:
            mi[com] = mi.get(com, 0.0) + edge.get(weight, 1.0)
    
    ed_measure = 0
    for com, mi_num in mi.items():
        ed_measure += mi_num
    
    return float(ed_measure / len(com2nodes[ci]))


def average_degree(graph, node2com, com2nodes, ci, weight='weight'):
    """
    Formula:
        return_value = (1.0 / ||C_i||) * SIGMA{ 2 * mi / ni }, whereas:
        ||C_i|| = ni : number of nodes in community ci
        mi : number of links in community ci

    Arguments:
        graph {[nx.Graph]} -- [network to be processed]
        node2com {[dict]} -- [keys are nodes, values are corresponding community index]
        com2nodes {[dict]} -- [keys are community indexed, values are nodes inside each community]
        ci {[int]} -- [index of the intended community]
    
    Keyword Arguments:
        weight {str} -- ['weight' or 'prob'] (default: {'weight'})
    
    Returns:
        [float] -- [average degree of community ci of the graph based on given partition]
    """
    mi = dict()
    for edge in graph.edges(data=True):
        com = node2com[edge[0]]
        if node2com[edge[1]] == com:
            mi[com] = mi.get(com, 0.0) + edge.get(weight, 1.0)
    
    ad_measure = 0
    for com, nodes in com2nodes.items():
        ad_measure += 2.0 * mi[com] / len(nodes)
    
    return float(ad_measure / len(com2nodes[ci]))


def FOMD(graph, node2com, com2nodes, ci, weight='weight'):
    """
    Formula:
        return_value = (1.0 / ||C_i||) * SIGMA{ (||A||) / ni }, whereas:
        ||C_i|| = ni : number of nodes in community ci
        A : {u,v in C_i , ||{u,v}|| > median_mi}, whereas:
        median_mi : median of degrees of nodes of the graph

    Arguments:
        graph {[nx.Graph]} -- [network to be processed]
        node2com {[dict]} -- [keys are nodes, values are corresponding community index]
        com2nodes {[dict]} -- [keys are community indexed, values are nodes inside each community]
        ci {[int]} -- [index of the intended community]
    
    Keyword Arguments:
        weight {str} -- ['weight' or 'prob'] (default: {'weight'})
    
    Returns:
        [float] -- [Fraction over median degree of community ci of the graph based on given partition]
    """
    degrees = dict(graph.degree(weight=weight))
    degree_values = list(degrees.values())
    median_mi = stat.median(degrees)

    com2edges = dict()
    for edge in graph.edges(data=True):
        com = node2com[edge[0]]
        if node2com[edge[1]] == com:
            com2edges[com] = com2edges.get(com, 0.0) + edge.get(weight, 1.0)

    fomd = 0
    for com, nodes in com2nodes.items():
        size_of_set = com2edges[com]
        if size_of_set > median_mi:
            fomd += float(size_of_set / len(nodes))

    return float(fomd / len(com2nodes[ci]))


def TPR(graph, node2com, com2nodes, ci):
    """ the fraction of nodes in C_i that belongs to a triad
    
    Arguments:
        graph {[nx.Graph]} -- [network to be processed]
        node2com {[dict]} -- [keys are nodes, values are corresponding community index]
        com2nodes {[dict]} -- [keys are community indexed, values are nodes inside each community]
        ci {[int]} -- [index of the intended community]
        
    Returns:
        [float] -- [triangle participation ratio]
    """
    triad_nodes = dict()
    for n1, n2, edge_info in graph.edges(data=True):
        com = node2com[n1]
        if node2com[n2] == com:
            n1_neighbors = set(graph.neighbors(n1))
            n2_neighbors = set(graph.neighbors(n2))
            mutual = list(n1_neighbors.intersection(n2_neighbors))
            for n3 in mutual:
                if node2com[n1] == node2com[n3]:
                    if graph.has_edge(n1, n3) and graph.has_edge(n2, n3):
                        triad_nodes.setdefault(com, set()).update([n1, n2, n3])
    tpr = 0.0
    for com, nodes in triad_nodes.items():
        tpr += (len(nodes) / com2nodes(com))

    return float(tpr / len(com2nodes[ci]))


def expansion(graph, node2com, com2nodes, ci, weight='weight'):
    """
    Formula:
        return_value = 1 / ||C_i|| * SIGMA { b_i / n_i }, whereas:
        b_i : edges that have one endpoint in the community c_i
        n_i : number of nodes in the community c_i
    
    Arguments:
        graph {[nx.Graph]} -- [network to be processed]
        node2com {[dict]} -- [keys are nodes, values are corresponding community index]
        com2nodes {[dict]} -- [keys are community indexed, values are nodes inside each community]
        ci {[int]} -- [index of the intended community]
    
    Keyword Arguments:
        weight {str} -- ['weight' or 'prob'] (default: {'weight'})
    
    Returns:
        [float] -- [number of edges per node that points outside the community]
    """
    bi = dict()
    for n1, n2, edge_info in graph.edges(data=True):
        comm1 = node2com[n1]
        comm2 = node2com[n2]
        if comm1 != comm2:
            bi[comm1] = bi.get(comm1, 0.0) + edge_info.get(weight, 1.0)
            bi[comm2] = bi.get(comm2, 0.0) + edge_info.get(weight, 1.0)
    
    e = 0.0
    for com, nodes in com2nodes.items():
        e += float(bi[com] / len(nodes))
    
    return float(e / len(com2nodes[ci]))


def cut_ratio(graph, node2com, com2nodes, ci, weight='weight'):
    """
    Formula:
        return_value = 1 / ||C_i|| * SIGMA{ b_i / (n_i * (n - n_i)) }, whereas:
        b_i : edges that have one endpoint in the community c_i
        n_i : number of nodes in the community c_i
        n : number of nodes in the network

    Arguments:
        graph {[nx.Graph]} -- [network to be processed]
        node2com {[dict]} -- [keys are nodes, values are corresponding community index]
        com2nodes {[dict]} -- [keys are community indexed, values are nodes inside each community]
        ci {[int]} -- [index of the intended community]
    
    Keyword Arguments:
        weight {str} -- ['weight' or 'prob'] (default: {'weight'})
    
    Returns:
        [float] -- [the fraction of existing edges leaving the community ci]
    """
    bi = dict()
    for n1, n2, edge_info in graph.edges(data=True):
        comm1 = node2com[n1]
        comm2 = node2com[n2]
        if comm1 != comm2:
            bi[comm1] = bi.get(comm1, 0.0) + edge_info.get(weight, 1.0)
            bi[comm2] = bi.get(comm2, 0.0) + edge_info.get(weight, 1.0)
    
    e = 0.0
    n = graph.number_of_nodes()
    for com, nodes in com2nodes.items():
        ni = len(nodes)
        e += float(bi[com] / (ni * (n - ni)))
    
    return float(e / len(com2nodes[ci]))


def conductance(graph, node2com, com2nodes, ci, weight='weight'):
    """
    Formula:
        return_value = 1 / ||C_i|| * SIGMA{ b_i / (2 * m_i + b_i) }, whereas:
        b_i : edges that have one endpoint in the community c_i
        m_i : number of links in the community c_i

    Arguments:
        graph {[nx.Graph]} -- [network to be processed]
        node2com {[dict]} -- [keys are nodes, values are corresponding community index]
        com2nodes {[dict]} -- [keys are community indexed, values are nodes inside each community]
        ci {[int]} -- [index of the intended community]
    
    Keyword Arguments:
        weight {str} -- ['weight' or 'prob'] (default: {'weight'})
    
    Returns:
        [float] -- [the fraction of total edges volume that points outside community ci]
    """
    bi = dict()
    mi = dict()
    for n1, n2, edge_info in graph.edges(data=True):
        comm1 = node2com[n1]
        comm2 = node2com[n2]
        if comm1 == comm2:
            mi[comm1] = mi.get(comm1, 0.0) + edge_info.get(weight, 1.0)
            mi[comm2] = mi.get(comm2, 0.0) + edge_info.get(weight, 1.0)
        else:
            bi[comm1] = bi.get(comm1, 0.0) + edge_info.get(weight, 1.0)
            bi[comm2] = bi.get(comm2, 0.0) + edge_info.get(weight, 1.0)

    c = 0.0
    for com, nodes in com2nodes.items():
        c += float(bi[com] / (2 * mi[com] + bi[com]))

    return float(c / len(com2nodes[ci]))


def normalized_cut(graph, node2com, com2nodes, ci, weight='weight'):
    """
    
    Arguments:
        graph {[nx.Graph]} -- [network to be processed]
        node2com {[dict]} -- [keys are nodes, values are corresponding community index]
        com2nodes {[dict]} -- [keys are community indexed, values are nodes inside each community]
        ci {[int]} -- [index of the intended community]
    
    Keyword Arguments:
        weight {str} -- ['weight' or 'prob'] (default: {'weight'})
    
    Returns:
        [float] -- [normalized cut]
    """
    bi = dict()
    mi = dict()
    for n1, n2, edge_info in graph.edges(data=True):
        comm1 = node2com[n1]
        comm2 = node2com[n2]
        if comm1 == comm2:
            mi[comm1] = mi.get(comm1, 0.0) + edge_info.get(weight, 1.0)
            mi[comm2] = mi.get(comm2, 0.0) + edge_info.get(weight, 1.0)
        else:
            bi[comm1] = bi.get(comm1, 0.0) + edge_info.get(weight, 1.0)
            bi[comm2] = bi.get(comm2, 0.0) + edge_info.get(weight, 1.0)

    nc = 0.0
    m = graph.size(weight=weight)
    for com, nodes in com2nodes.items():
        nc += bi[com] * (1.0 / (2 * mi[com] + bi[com]) + 1.0 / (2 * (m - mi[com]) + bi[com]))

    return float(nc / len(com2nodes[ci]))


def MODF(graph, node2com, com2nodes, ci, weight='weight'):
    """
    
    Arguments:
        graph {[nx.Graph]} -- [network to be processed]
        node2com {[dict]} -- [keys are nodes, values are corresponding community index]
        com2nodes {[dict]} -- [keys are community indexed, values are nodes inside each community]
        ci {[int]} -- [index of the intended community]
    
    Keyword Arguments:
        weight {str} -- ['weight' or 'prob'] (default: {'weight'})
    
    Returns:
        [float] -- [maximum out degree fraction]
    """
    out_neighs = dict()
    for n1, n2, edge_info in graph.edges(data=True):
        comm1 = node2com[n1]
        comm2 = node2com[n2]
        if comm1 != comm2:
            mi[n1] = out_neighs.get(comm1, 0.0) + edge_info.get(weight, 1.0)
            mi[n2] = out_neighs.get(comm2, 0.0) + edge_info.get(weight, 1.0)

    modf_dict = dict()
    for com, nodes in com2nodes.items():
        max_ratio = 0.0
        for node in nodes:
            ratio = float(out_neighs.get(node, 0.0) / graph.degree[node])
            if ratio > max_ratio:
                max_ratio = ratio
        modf_dict[com] = max_ratio

    modf = 0.0
    for com, max_ratio in modf_dict.items():
        modf += max_ratio

    return float(modf / len(com2nodes[ci]))


def AODF(graph, node2com, com2nodes, ci, weight='weight'):
    """
    
    Arguments:
        graph {[nx.Graph]} -- [network to be processed]
        node2com {[dict]} -- [keys are nodes, values are corresponding community index]
        com2nodes {[dict]} -- [keys are community indexed, values are nodes inside each community]
        ci {[int]} -- [index of the intended community]
    
    Keyword Arguments:
        weight {str} -- ['weight' or 'prob'] (default: {'weight'})
    
    Returns:
        [float] -- [average out degree fraction]
    """
    out_neighs = dict()
    for n1, n2, edge_info in graph.edges(data=True):
        comm1 = node2com[n1]
        comm2 = node2com[n2]
        if comm1 != comm2:
            mi[n1] = out_neighs.get(comm1, 0.0) + edge_info.get(weight, 1.0)
            mi[n2] = out_neighs.get(comm2, 0.0) + edge_info.get(weight, 1.0)

    modf_dict = dict()
    for com, nodes in com2nodes.items():
        for node in nodes:
            ratio = float(out_neighs.get(node, 0.0) / graph.degree[node])
            modf_dict[com] = modf_dict.get(com, 0.0) + ratio

    modf = 0.0
    for com, ratio in modf_dict.items():
        modf += float(ratio / len(com2nodes[com]))

    return float(modf / len(com2nodes[ci]))


def FODF(graph, node2com, com2nodes, ci, weight='weight'):
    """
    
    Arguments:
        graph {[nx.Graph]} -- [network to be processed]
        node2com {[dict]} -- [keys are nodes, values are corresponding community index]
        com2nodes {[dict]} -- [keys are community indexed, values are nodes inside each community]
        ci {[int]} -- [index of the intended community]
    
    Keyword Arguments:
        weight {str} -- ['weight' or 'prob'] (default: {'weight'})
    
    Returns:
        [float] -- [flake out degree fraction]
    """
    node_dict = dict()
    for node in graph.nodes():
        comm = node2com[node]
        deg = graph.degree[node]
        neighbors = list(graph.neighbors(node))
        for neigh in neighbors:
            if node2com[neigh] == comm:
                node_dict[node] = node_dict.get(node, 0) + graph[node][neigh].get(weight, 1.0)
        node_dict[node] = node_dict.get(node, 0.0) - 0.5 * deg

    fodf = 0.0
    for com, nodes in com2nodes.items():
        fodf += len([x for x in nodes if node_dict[node] < 0.0]) / len(nodes)

    return float(fodf / len(com2nodes[ci]))


def modularity(graph, com2nodes, la=1.0, weight='weight'):
    """
    
    Arguments:
        graph {[nx.Graph]} -- [network to be processed]
        com2nodes {[dict]} -- [keys are community indexed, values are nodes inside each community]
    
    Keyword Arguments:
        la {float} -- [lambda parameter in the formula] (default: {1.0})
        weight {str} -- ['weight' or 'prob'] (default: {'weight'})
    
    Returns:
        [float] -- [quality of the graph based on modularity introduced by Newman]
    """
    q = 0.0
    m = graph.size(weight=weight)
    degrees = dict(graph.degree(weight=weight))
    for com, nodes in com2nodes.items():
        for (n1, n2) in [(x, y) for x in nodes for y in nodes]:
            w = 0.0
            if graph.has_edge(n1, n2):
                w = graph[n1][n2].get(weight, 1.0)
            if n1 == n2:
                w *= 2
            q = q + w - la * degrees[n1] * degrees[n2] / (2 * m)

    return q / (2 * m)


def NMI(node2com_true_filename, node2com_pred):
    labels_true_dict = dict()
    with open(node2com_true_filename, 'r') as file:
        lines = file.readlines()
        for line in lines:
            line = line.split()
            labels_true_dict[int(line[0])] = int(line[1])

    labels_true = list()
    labels_pred = list()

    for key in labels_true_dict.keys() & node2com_pred.keys():
        labels_true.append(labels_true_dict[key])
        labels_pred.append(node2com_pred[key])

    return normalized_mutual_info_score(labels_true, labels_pred, average_method='arithmetic')