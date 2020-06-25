import utils
import time
import networkx as nx
import random
import os.path
from copy import deepcopy
from measures import modularity
from measures import NMI


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
    if Tin > Tout:
        return Tin * (Tin - Tout)
    return 0.0


def find_best_next_node(graph, community, shell_set, Tin_prev, Tout_prev):
    dict_of_Tin = dict()
    dict_of_Tex = dict()
    dict_of_T = dict()

    for node in shell_set:
        T_in_node = calc_T_in(graph, node, community, Tin_prev)
        T_out_node = calc_T_out(graph, node, community, Tout_prev)
        dict_of_Tin[node] = T_in_node
        dict_of_Tex[node] = T_out_node
        dict_of_T[node] = calc_T(T_in_node, T_out_node)

    best_node = shell_set[0]
    for i in range(1, len(shell_set)):
        if dict_of_T[shell_set[i]] > dict_of_T[best_node]:
            best_node = shell_set[i]
        elif dict_of_T[shell_set[i]] == dict_of_T[best_node]:
            if dict_of_Tex[shell_set[i]] < dict_of_Tex[best_node]:
                best_node = shell_set[i]

    return best_node, dict_of_Tin[best_node], dict_of_Tex[best_node]


def find_node_highest_degree(graph, communities):
    nodes = list(graph.nodes())

    for comm in communities:
        for node in comm:
            if node in nodes:
                nodes.remove(node)

    highest_deg_node = nodes[0]
    for i in range(1, len(nodes)):
        if graph.degree[nodes[i]] > graph.degree[highest_deg_node]:
            highest_deg_node = nodes[i]

    return highest_deg_node


def community_search(graph):    
    communities = list()
    nodes_discovered = 0
    while nodes_discovered < graph.number_of_nodes():
        initial_node = find_node_highest_degree(graph, communities)
        community = [initial_node]
        nodes_discovered += 1
        
        shell_set = list(graph.neighbors(initial_node))
        for com in communities:
            for node in com:
                if node in shell_set:
                    shell_set.remove(node)

        Tin = calc_T_in(graph, initial_node, [], 0.0)
        Tout = calc_T_out(graph, initial_node, [], 0.0)

        while len(community) < graph.number_of_nodes():
            best_node, new_Tin, new_Tout = find_best_next_node(graph, community, shell_set, Tin, Tout)
            T = calc_T(new_Tin, new_Tout)
            if T >= calc_T(Tin, Tout):
                Tin, Tout = new_Tin, new_Tout
                new_neighbors = list(set(graph.neighbors(best_node)) - set(community))
                if len(communities) > 0:
                    for comm in communities:
                        for node in comm:
                            if node in new_neighbors:
                                new_neighbors.remove(node)

                community.append(best_node)
                nodes_discovered += 1

                shell_set.extend(new_neighbors)
                shell_set = list(set(shell_set))
                shell_set.remove(best_node)

                if not shell_set:
                    break
            
            else:
                break

        communities.append(community)

    return sorted(communities)



def main():
    start_time = time.time()

    args = utils.create_argument_parser()
    graph = utils.load_graph(args.dataset, args.w)
    graph_copy = deepcopy(graph)
    
    communities = community_search(graph)
    com_dict = {}

    for i in range(len(communities)):
        com_dict[i] = communities[i]

    com_dict2 = {}
    for k, v in com_dict.items():
        for node in v:
            com_dict2[node] = k

    utils.print_comm_info_to_display(com_dict)

    print('modularity_value =', modularity(graph, com_dict))            
    print('NMI =', NMI(args.output, com_dict2))

    finish_time = time.time()
    print('\nDone in %.4f seconds.' %(finish_time - start_time))


if __name__ == "__main__":
    main()