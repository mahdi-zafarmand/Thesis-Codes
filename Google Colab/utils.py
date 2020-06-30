import argparse
import os.path
import networkx as nx
from collections import defaultdict
import measures


def create_argument_parser():
	"""This function reads the command line argument and parse them
	
	Returns:
		[parser] -- [the parsed arguments in received in the command line]
	"""
	parser = argparse.ArgumentParser(description='The SIWO algorithm', prog="SIWO")
	parser.add_argument('address1', help='Path to the network1 files')
	parser.add_argument('address2', help='Path to the network1 files')
	parser.add_argument('-w', help='weighted graphs', action="store_true")
	parser.add_argument('-r', help='randomize order of nodes', action="store_true")
	return parser.parse_args()


def load_graph(path, weighted, self_loop=False):
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
			if self_loop:
				graph.add_edge(v1, v2, weight = w)                
			elif v1 != v2:
				graph.add_edge(v1, v2, weight = w)

	# Build the largest connected component if the whole graph is not connected.
	if not nx.is_connected(graph):
		print('graph is not originally connected.')
		nodes_set = max(nx.connected_components(graph), key=len)
		connected_subgraph = graph.subgraph(nodes_set)
		graph = nx.Graph(graph) # creating a graph out of the extracted subgraph to avoid the issue with frozen graphs.
		print('The graph has', graph.number_of_nodes(), ' nodes and is completely loaded!')

	return graph


def extract_communities(node2com):
	"""converts node2com to com2nodes

	Args:
		node2com ([dict]): [a map of nodes to their corresponding community indices]

	Returns:
		[dict]: [a map of community indices to their corresponding nodes]
	"""
	com2nodes = defaultdict(list)
	for node, com in node2com.items():
		com2nodes[com].append(node)

	return dict(com2nodes)


def find_node_highest_degree(graph, communities=[]):
	"""finds the node with highest degree among all node of the network that are not discovered yet.

	Args:
		graph ([nx.Graph]): [the given network]
		communities (list, optional): [all discovered communities]. Defaults to [].

	Returns:
		[int]: [a node with the highest degree]
	"""
	nodes = list(graph.nodes())

	for community in communities:
		for node in community:
			if node in nodes:
				nodes.remove(node)

	highest_deg_node = nodes[0]
	for i in range(1, len(nodes)):
		if graph.degree[nodes[i]] > graph.degree[highest_deg_node]:
			highest_deg_node = nodes[i]

	return highest_deg_node


def find_all_new_candidates(graph, community):
	"""finds all nodes that are directly connected to at least one member of the community.

	Args:
		graph ([nx.Graph]): [the given network]
		community ([list]): [contains nodes that are already discovered to be in the community]

	Returns:
		[list]: [all direct neighbors of nodes in the community which are not in the community]
	"""
	candidates = set()
	for node in community:
		for neigh in graph.neighbors(node):
			candidates.add(neigh)

	return list(candidates - set(community))


def clear_former_nodes(candidates, communities):
	"""clear discoevered nodes from the list of candidates.

	Args:
		candidates ([list]): [list of candidates]
		communities ([list]): [list of all already discovered communities]

	Returns:
		[list]: [list of proper candidates after clearing up extras]
	"""
	for community in communities:
		for node in community:
			if node in candidates:
				candidates.remove(node)
	return candidates


def read_ground_truth(file_address):
	"""reads the ground-truth file and creates a dict in which keys are community indices and values are their corresponding nodes.

	Args:
		file_address ([str]): [filename of the ground-truth information of communities]

	Returns:
		[dict]: [keys: community indices and values: their corresponding nodes]
	"""
	com2nodes = defaultdict(list)
	with open(file_address, 'r') as file:
		lines = file.readlines()
		for line in lines:
			line = line.split()
			com2nodes[int(line[1])].append(int(line[0]))

	return com2nodes


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

	print('\nSize of communities are:', sizes)
	print('There are', len(sizes), 'communities detected.')


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


def remove_self_loops(graph):
	"""removes the edges that connects an edge in the network to itself

	Args:
		graph ([nx.Graph]): [the given network]

	Returns:
		[nx.Graph]: [the same graph without any self loops]
	"""
	for node in graph.nodes():
		if graph.has_edge(node, node):
			graph.remove_edge(node, node)
	return graph


def report_performance(graph, partition, ground_truth_file_address):
	"""prints some performace measure of the experiment, inclusing Qmodularity and NMI score of the partition.

	Args:
		graph ([nx.Graph]): [the given network]
		partition ([list]): [list of all discovered communities]
		ground_truth_file_address ([str]): [filename of the ground-truth information of communities]
	"""
	com2nodes = dict()
	for i in range(len(partition)):
		com2nodes[i] = partition[i]

	node2com = dict()
	for com_index, nodes in com2nodes.items():
		for node in nodes:
			node2com[node] = com_index

	report_str = 'Modularity = ' +str(measures.modularity(graph, com2nodes))[:5] + '\t'
	report_str += 'NMI = ' +str(measures.NMI(ground_truth_file_address, node2com))[:5]
	return report_str


def update_performance_info(node, performance_info, community_pred, ground_truth_com2nodes):
	"""updates the performance measures for a node based on the ground-truth community information and detected community.

	Args:
		node ([int]): [a node of the network]
		performance_info ([dict]): [contains precision, recall, and f1-score for each node]
		community_pred ([list]): [list of all nodes in the same community as of $node]
		ground_truth_com2nodes ([dict]): [keys: nodes of the network, values: nodes of the corresponding community index]
	"""
	community_true = list()
	for com_index, nodes in ground_truth_com2nodes.items():
		if node in nodes:
			community_true = list(nodes)
			break
	
	num_common_elements = len(set(community_true).intersection(set(community_pred)))

	try:
		recall = num_common_elements / len(community_true)
	except:
		recall = 0.0
	try:
		precision = num_common_elements / len(community_pred)
	except:
		precision = 0.0
	try:
		f1_score = 2 * (recall * precision) / (recall + precision)
	except:
		f1_score = 0.0

	performance_info[node]['recall'] = recall
	performance_info[node]['precision'] = precision
	performance_info[node]['f1_score'] = f1_score


def amend_community_size_one(graph, partition, size_one_coms):
	"""merges given communities of size one into another community that suits the best.

	Args:
		graph ([nx.Graph]): [the given network]
		partition ([list]): [list of all discovered communities with more than two nodes]
		size_one_coms ([list]): [list of network's all communities of size one]

	Returns:
		[list]: [list of all communities after amending the communities of size one]
	"""
	for com in size_one_coms:
		neighbors = set(graph.neighbors(com[0]))
		# define dict of strength for comparison
		strength = dict()
		for neigh in neighbors:
			for i in range(len(partition)):
				if neigh in partition[i] and graph.has_edge(com[0], neigh):
					strength[i] = strength.get(i, 0.0) + graph[com[0]][neigh].get('weight', 0.0)
					break

		# find best community based on highest strength
		best_com_index = list(strength.keys())[0]
		for i in strength:
			if strength[i] > strength[best_com_index]:
				best_com_index = i
		partition[best_com_index] = sorted(partition[best_com_index] + com)

	return partition


def amend_community_size_two(graph, partition, size_two_coms):
	"""merges given communities of size two into another community that suits the best.

	Args:
		graph ([nx.Graph]): [the given network]
		partition ([list]): [list of all discovered communities with more than two nodes]
		size_two_coms ([list]): [list of network's all communities of size two]

	Returns:
		[list]: [list of all communities after amending the communities of size two]
	"""
	for com in size_two_coms:
		neighbors = set(graph.neighbors(com[0]))
		neighbors.update(graph.neighbors(com[1]))
		# define dict of strength for comparison
		strength = dict()
		for neigh in neighbors:
			for i in range(len(partition)):
				if neigh in partition[i]:
					if graph.has_edge(com[0], neigh):
						strength[i] = strength.get(i, 0.0) + graph[com[0]][neigh].get('weight', 0.0)
					if graph.has_edge(com[1], neigh):
						strength[i] = strength.get(i, 0.0) + graph[com[1]][neigh].get('weight', 0.0)
					break

	# find best community based on highest strength
		best_com_index = list(strength.keys())[0]
		for i in strength:
			if strength[i] > strength[best_com_index]:
				best_com_index = i
		partition[best_com_index] = sorted(partition[best_com_index] + com)

	return partition


def amend_partition(graph, partition):
	"""amends the detected partition by merging communities of size 1 or 2 to the other communities in the partition.

	Args:
		graph ([nx.Graph]): [the given network]
		partition ([list]): [list of all detected communities]

	Returns:
		[list]: [list of all communities after amendment]
	"""
	size_one_coms = [x for x in partition if len(x) == 1]
	size_two_coms = [x for x in partition if len(x) == 2]

	# remove communities of size 1 or 2 from the partition
	if size_one_coms or size_two_coms:
		i = 0
		while i < len(partition):
			if len(partition[i]) < 3:
				partition.pop(i)
				i -= 1
			i += 1

	if size_one_coms != []:
		partition = amend_community_size_one(graph, partition, size_one_coms)

	if size_two_coms != []:
		partition = amend_community_size_two(graph, partition, size_two_coms)

	return partition