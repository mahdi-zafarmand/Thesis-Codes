import utils
import time
import networkx as nx
import random
import os.path


def find_initial_node(graph, node):
    neighbors = list(graph.neighbors(node))
    neighbors.append(node)
    
    deg = 0.0
    candidates = list()
    for neigh in neighbors:
        deg_compare = graph.degree[neigh]
        if deg_compare == deg:
            candidates.append(neigh)
        elif deg_compare > deg:
            deg = deg_compare
            candidates.clear()
            candidates.append(neigh)

    return random.choice(candidates)


def calc_T_in(graph, node, community, Tin_prev):
    increment = 0.0
    neighbors = list(graph.neighbors(node))

    for neigh1 in neighbors:
        for neigh2 in neighbors:
            if neigh1 != neigh2 and graph.has_edge(neigh1, neigh2):
                if neigh1 in community and neigh2 in community:
                    increment += 1

    return Tin_prev + 0.5 * increment


def calc_T_out(graph, node, community, Tout_prev):
    increment = 0.0
    decrement = 0.0
    neighbors = list(graph.neighbors(node))

    for neigh1 in neighbors:
        for neigh2 in neighbors:
            if neigh1 != neigh2 and graph.has_edge(neigh1, neigh2):
                if not neigh1 in community and not neigh2 in community:
                    increment += 1
                elif neigh1 in community and not neigh2 in community:
                    decrement += 1

    return Tout_prev + 0.5 * increment - decrement


def calc_T(Tin, Tout):
    if Tin> Tout:
        return Tin * (Tin - Tout)
    return 0


def run_one_iter(graph, community, shell_set, Tin_prev, Tout_prev):
    dict_of_Tin = dict()
    dict_of_Tex = dict()
    dict_of_T = dict()

    for node in shell_set:
        T_in_node = calc_T_in(graph, node, community, Tin_prev)
        T_out_node = calc_T_out(graph, node, community, Tout_prev)
        dict_of_Tin[node] = T_in_node
        dict_of_Tex[node] = T_out_node
        dict_of_T[node] = calc_T(T_in_node, T_out_node)

    best_value = 0.0
    best_nodes = list()
    for node, value in dict_of_T.items():
        if value == best_value:
            best_nodes.append(node)
        elif value > best_value:
            best_value = value
            best_nodes.clear()
            best_nodes.append(node)
    
    lowest_tx = float('inf')
    best_nodes_tie_break = list()
    for node, value in dict_of_Tex.items():
        if value == lowest_tx:
            best_nodes_tie_break.append(node)
        elif value < lowest_tx:
            lowest_tx = value
            best_nodes_tie_break.clear()
            best_nodes_tie_break.append(node)

    best_node = random.choice(best_nodes_tie_break)
    return best_node, dict_of_Tin[best_node], dict_of_Tex[best_node]


def find_next_comm(graph):
    nodes = list(graph.nodes())
    random_initial_node = random.choice(nodes)
    initial_node = find_initial_node(graph, random_initial_node)
    # print('random node =', random_initial_node, ' starting node =', initial_node)

    community = [initial_node]
    shell_set = list(graph.neighbors(initial_node))
    Tin = calc_T_in(graph, initial_node, [], 0.0)
    Tout = calc_T_out(graph, initial_node, [], 0.0)

    while len(community) < len(nodes):
        best_node, new_Tin, new_Tout = run_one_iter(graph, community, shell_set, Tin, Tout)
        T = calc_T(new_Tin, new_Tout)
        if T >= calc_T(Tin, Tout):
            Tin, Tout = new_Tin, new_Tout
            new_neighbors = list(set(graph.neighbors(best_node)) - set(community))    
            community.append(best_node)
            shell_set.extend(new_neighbors)
            shell_set = list(set(shell_set))
            shell_set.remove(best_node)

            if not shell_set:
                break
        
        else:
            break
        
    return community


def find_comms(graph):
    partition = list()
    while graph.number_of_nodes() > 0:
        comm = find_next_comm(graph)
        comm.sort()
        partition.append(comm)
        print(len(partition), ' : ', comm, ' len =', len(comm))
        for node in comm:
            if node in graph.nodes():
                graph.remove_node(node)

    return partition


def load_graph(path, weighted, report=True):
    """This function reads the edge list then make the graph.
    
    Arguments:
        path {[str]} -- [path of the file containing the edge list information]
        weighted {[bool]} -- [shows if a graph should be considered as weighted or not]
        
    Returns:
        [nx.Graph] -- [a graph made based on the edge list information]
    """
    graph = nx.Graph()
    # Read the graph
    if(not os.path.isfile(path)):
        print("Error: file " + path + "not found!")
        exit(-1)
    with open(path) as f:
        delimiter = '\t'
        for line in f.readlines():
            line = line.split(delimiter)
            v2 = int(line[0])
            v1 = int(line[1])
            p = float(line[2])
            graph.add_node(v1, size = 1)
            graph.add_node(v2, size = 1)
            if v1 != v2:
                graph.add_edge(v1, v2, prob = p)

    # Build the largest connected component if the whole graph is not connected.
    if not nx.is_connected(graph):
        print('graph is not originally connected.')
        nodes_set = max(nx.connected_components(graph), key=len)
        connected_subgraph = graph.subgraph(nodes_set)
        graph = nx.Graph(graph) # creating a graph out of the extracted subgraph to avoid the issue with frozen graphs.
    
    if report:
        print('graph has', len(graph.nodes()), ' nodes.')
        print('Graph is completely loaded!')

    return graph


def main():
    start_time = time.time()

    args = utils.create_argument_parser()
    graph = utils.load_graph(args.dataset, args.w)

    partition = find_comms(graph)

    finish_time = time.time()
    print('\nDone in %.4f seconds.' %(finish_time - start_time))


if __name__ == "__main__":
    main()