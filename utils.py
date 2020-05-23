import argparse
import os.path
import networkx as nx
from collections import defaultdict

def create_argument_parser():
    """This function reads the command line argument and parse them
    
    Returns:
        [parser] -- [the parsed arguments in received in the command line]
    """
    parser = argparse.ArgumentParser(description='The SIWO algorithm', prog="SIWO")
    parser.add_argument('dataset', help='Path to the network files')
    parser.add_argument('output', help='output name')
    parser.add_argument('-w', help='weighted graphs', action="store_true")
    parser.add_argument('-r', help='randomize order of nodes', action="store_true")
    return parser.parse_args()


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
            if line.find(delimiter) < 0:
                delimiter = ' '
            w = 1.0
            line = line.split(delimiter)
            if weighted:
                w = float(line[2])
            v2 = int(line[0])
            v1 = int(line[1])
            graph.add_node(v1, size = 1)
            graph.add_node(v2, size = 1)
            if v1 != v2:
                graph.add_edge(v1, v2, weight = w)

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


def extract_communities(node2com):
	com2nodes = defaultdict(list)
	for node, com in node2com.items():
		com2nodes[com].append(node)

	return dict(com2nodes)


def print_comm_info_to_display(communities, size_info=True, lone_comms = True):
    """This function displays information about the partitions.
    
    Arguments:
        partition {[dict]} -- [contains information about the communities]
    
    Keyword Arguments:
        lone_comms {bool} -- [if you want to show the lone communities or not] (default: {True})
    """

    for comm_index in communities:
        comm_size = len(communities[comm_index])
        if not lone_comms and comm_size == 1:
            continue
        print(comm_index, '-> size =', comm_size, '-> members =', sorted(communities[comm_index]))

    sizes = []
    for comm_index in communities:
        sizes.append(len(communities[comm_index]))
    sizes.sort()

    print('\nsize of communities are:', sizes)
    print('there are', len(sizes), 'communities detected.')


def write_comm_info_to_file(output, partition):
    """This function writes the results to a file.
    
    Arguments:
        output {[str]} -- [the path to the output file]
        partition {[dict]} -- [dictionary contains information about the communities.]
    """
    f = open("results/" + output, 'w')
    for elem, part in sorted(partition.items()):
        out = str(elem) + " " + str(part)
        f.write(out + "\n")
    f.close()


