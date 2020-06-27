import networkx as nx
import numpy as np
import random
import measures
import time


def generate_the_network(n, tau1, tau2, mu, min_degree, max_degree):
	"""generates the network with the given parameter using the LFR benchmark 
	method from the networkx package. This function tries 10 times to generate
	the network, if it cannot it returns None.

	Args:
		n ([int]): [number of nodes in the network]
		tau1 ([float]): [power-law exponent for degree distribution]
		tau2 ([float]): [power-law exponent for community size distribution]
		mu ([float]): [fraction of intra-community edges to each node]
		min_degree ([int]): [minimum degree of a node in the network]
		max_degree ([int]): [maximum degree of a node in the network]

	Returns:
		[nx.Graph]: [synthetic network]
	"""
	num_of_tries = 10
	successfully_built = False
	while num_of_tries > 0 and successfully_built == False:
		# print('i =', i)   # for debug purposes only!
		try:
			G = nx.LFR_benchmark_graph(n, tau1, tau2, mu, min_degree=min_degree, max_degree=max_degree)
			successfully_built = True
		except Exception as e:
			pass
		num_of_tries -= 1

	if successfully_built == False:
		return None

	return G


def extract_communities(G):
	"""extracts the groun-truth communities from the generated graph G.

	Args:
		G ([nx.Graph]): [the corresponding graph of the desired network]

	Returns:
		[tuple]: [contains two sets of maps, one from nodes to their community indices
		and the other from community indices to their nodes]
	"""
	communities = {frozenset(G.nodes[v]['community']) for v in G}
	com2nodes = dict()
	node2com = dict()
	n = 0
	for nodes_in_community in communities:
		com2nodes[n] = list(nodes_in_community)
		for node in nodes_in_community:
			node2com[node] = n
		n += 1
	return node2com, com2nodes


def print_graph_info(G, com2nodes, decimal_precision = 3, comm_size_info = False):
	"""prints some information about the network, including the number of nodes, number of edges,
	the average degree, the modularity value, number and sizes of the communities.

	Args:
		G ([nx.Graph]): [the given network]
		com2nodes ([dict]): [a map of community indices to their corresponding nodes]
		decimal_precision (int, optional): [determines the precision of shown float numbers]. Defaults to 3.
		comm_size_info (bool, optional): [determines if community size info should be printed]. Defaults to False.
	"""
	print('|V| =', G.number_of_nodes(), '\t|E| =', G.number_of_edges(), '\t|C| =', len(com2nodes), '\tAvg(degree) =', sum(G.degree[x] for x in G.nodes())/G.number_of_nodes())
	print('Modularity based on Newman formula =', round((measures.modularity(G, com2nodes)), decimal_precision))

	if comm_size_info == True:
		print('sizes of comms = [', end="")
		for i in range(len(com2nodes)):
			if i == len(com2nodes) - 1:
				print(len(com2nodes[i]), end="")
			else:
				print(len(com2nodes[i]), ', ', end="")
		print(']\n')


def write_to_file(G, node2com, com2nodes, edgelist_filename, community_filename):
	"""creates two files one containing a list of all edges and the other one containing a map of nodes
	to their ground-truth community index

	Args:
		G ([nx.Graph]): [the given network]
		node2com ([dict]): [a map of nodes to their corresponding community indices]
		com2nodes ([dict]): [a map of community indices to their corresponding nodes]
	"""
	num_nodes = G.number_of_nodes()
	nx.write_edgelist(G, edgelist_filename, delimiter='\t', data=False)

	with open(community_filename, 'w') as file:
		for node, com in sorted(node2com.items()):
			lineTowrite = str(node) + '\t' + str(com) + '\n'
			file.write(lineTowrite)
	

def run_generator(num_nodes, tau1, tau2, mu, min_degree):
	"""runs the network generator, extract ground-truth communities, print network information, and write to file.

	Args:
		num_nodes ([int]): [number of nodes in the network]
		tau1 ([float]): [power-law exponent for degree distribution]
		tau2 ([float]): [power-law exponent for community size distribution]
		mu ([float]): [fraction of intra-community edges to each node]
		min_degree ([int]): [minimum degree of a node in the network]
		max_degree ([int]): [maximum degree of a node in the network]
	"""
	start_time = time.time()

	graph_name = 'edgelist_' + str(num_nodes) + '_' + str(tau1) + '_' + str(tau2) + '_' + str(mu) + '_' + str(min_degree)
	comms_name = 'comslist_' + str(num_nodes) + '_' + str(tau1) + '_' + str(tau2) + '_' + str(mu) + '_' + str(min_degree)

	synthetic_network = generate_the_network(num_nodes, tau1, tau2, mu, min_degree, max_degree=num_nodes)
	if synthetic_network != None:
		print('Network generated successfully!')
		node2com, com2nodes = extract_communities(synthetic_network)
		graph_name += '_' + str(len(com2nodes)) + '.txt'
		comms_name += '_' + str(len(com2nodes)) + '.txt'

		print_graph_info(synthetic_network, com2nodes)
		write_to_file(synthetic_network, node2com, com2nodes, graph_name, comms_name)

	finish_time = time.time()
	print('-> Graph generated in %.4f seconds.' %(finish_time - start_time))
	return graph_name, comms_name

