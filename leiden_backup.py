import networkx as nx
from numpy.random import choice
from random import shuffle
from math import exp
from numpy import argmax
import utils
from copy import deepcopy
from itertools import product


# PREDETERMINED AND CONSTANT PARAMETERS OF THE METHOD
__GAMMA = 1.0					# the resolution parameter of the algorithm
__THETA = 0.01					# the exponential parameter of the random selection
__PREDEFINED_INF = 10 ** 20		# to avoid "max range error" in random selection
__IMPROV_THRESHOLD = 0.01		# to stop the iterations if no improvement is achieved
__ITER_THRESHOLD = 10			# to step the iterations based on number of iterations


def modularity(graph, partition):
	"""This function computes the modularity value based on the constant-potts model used in the Leiden paper.

	Arguments:
		graph {[nx.Graph]} -- [network that its single partitions will be found]
		partition {[list]} -- [contains lists that each one represents a community in the network]

	Returns:
		[float] -- [the modularity value]
	"""
	mod = 0.0
	for node1, node2, edge_info in graph.edges(data=True):
		comm_index1 = FindNodeInPartition(partition, node1)
		comm_index2 = FindNodeInPartition(partition, node2)
		if comm_index1 == comm_index2:
			mod += edge_info.get('weight', 1.0)


	m = max(1, graph.size(weight='weight'))
	for nodes in partition:
		size_comm = RecursiveSize(graph, nodes)
		mod -= (__GAMMA / (2.0 * m)) * (size_comm * (size_comm - 1) / 2.0)

	print()
	return mod


def SingletonPartition(graph):
	""" This function creates a new partition and puts every node of the given graph in a separate community.
	The main output is a list named "partition" and it has n communities which are lists themeselves.
	Here, n = |V|, if V is the set of nodes of the given graph. Every community has only one node in it.
	
	Arguments:
		graph {[nx.Graph]} -- [network that its single partitions will be found]
	
	Returns:
		[list] -- [contains lists with size one that each one represents a community in the network]
	"""
	partition = list()
	for node in graph.nodes():
		partition.append([node])
	return sorted(partition)


def FindNodeInPartition(all_communities, node):
	"""This auxiliary function finds the index of community in all_communities (the partition) that the given node belongs to.
	If (-1) is returned something is wrong!
	
	Arguments:
		all_communities {[list]} -- [contains lists that each one represents a community in the network]
		node {[int]} -- [node in the network]
	
	Returns:
		[int] -- [the index of desired community in all_communities]
	"""
	for comm_index in range(len(all_communities)):
		if node in all_communities[comm_index]:
			return comm_index

	print('SOMETHING IS WRONG IN **AggregateGraphHelper** FUNCTION.')
	return -1


def AggregateGraph(graph, partition, old_p):
	"""This function aggregates the former network to create the new one out of it, such that each node of the
	new graph represents a community in the former one. Since the size of the former communities are required 
	to compute the quality function, a node attribute is granted named "size" which represents how many nodes 
	of the given graph are inside a node of the new one. 
	
	Arguments:
		graph {[nx.Graph]} -- [network that will be aggregated]
		partition {[list]} -- [contains lists that each one represents a community in the network]
	
	Returns:
		[nx.Graph] -- [the aggregated network]
	"""
	# print(old_p)
	V_size_attr = dict()
	V_comm_info = dict()
	E_weight_attr = dict()

	# iterate over all communities to get their size
	for comm_index in range(len(partition)):
		V_size_attr[comm_index] = RecursiveSize(graph, partition[comm_index])
		V_comm_info[comm_index] = FindNodeInPartition(old_p, partition[comm_index][0])

	# iterate over all edges of the given graph to compute the new weights of the edges of the new aggregated graph
	for node1, node2, edge_info in graph.edges(data=True):
		comm_index1 = FindNodeInPartition(partition, node1)
		comm_index2 = FindNodeInPartition(partition, node2)
		if comm_index1 > comm_index2:
			comm_index1, comm_index2 = comm_index2, comm_index1
		if comm_index1 >= 0 and comm_index2 >= 0:
			E_weight_attr[(comm_index1, comm_index2)] = E_weight_attr.get((comm_index1, comm_index2), 0) + edge_info.get('weight', 1.0)
	
	# creating the new graph and assign nodes, edges and proper attributes to it.
	new_aggregated_graph = nx.Graph()
	new_aggregated_graph.add_nodes_from(V_size_attr)
	new_aggregated_graph.add_edges_from(E_weight_attr)
	nx.set_node_attributes(new_aggregated_graph, V_size_attr, name = 'size')
	nx.set_node_attributes(new_aggregated_graph, V_comm_info, name = 'comm')
	nx.set_edge_attributes(new_aggregated_graph, E_weight_attr, name = 'weight')

	return new_aggregated_graph


def flat(S):
	"""This function finds all nodes in the given list S or in its items (if their are lists too), 
	then puts all of the found items together in the returned parameter "flatten".
	
	Arguments:
		S {[list]} -- [may have nodes or list of nodes in it]
	
	Returns:
		[list] -- [the flattened list, only have nodes in it]
	"""
	S_copy = deepcopy(S)
	flatten = list()

	for item in S_copy:
		if type(item) == list:
			x = flat(item)
			flatten.extend(x)
		else:
			flatten.extend([item])
	
	return flatten


def RecursiveSize(graph, S):
	"""calculates the size of S based on the definittion in the paper, the number of nodes in S.
	
	Arguments:
		graph {[nx.Graph]} -- [network that is being investigated]
		S {[list]} -- [may have nodes or list of nodes inside it]
	
	Returns:
		[int] -- [the sum of sizes of recursively found nodes in S]
	"""
	s_flatten = flat(S)

	special_size_val = 0
	for node in s_flatten:
		special_size_val += graph[node].get('size', 1)

	return special_size_val


def RefinePartition(graph, partition):
	"""This function moves nodes in the already found communities and makes the refined partitions.
	
	Arguments:
		graph {[nx.Graph]} -- [the network that will be processed]
		partition {[list]} -- [contains lists that each one represents a community in the network]
	
	Returns:
		[list] -- [contains lists that each one represents a community in the network]
	"""
	partition_copy = deepcopy(partition)

	initial_modularity = modularity(graph, partition)

	# first, every node goes to its own community of size one.
	p_refined = SingletonPartition(graph)

	# the refinement is done for each found community separately, then the p_refined gets updated.
	for community_of_main_partition in partition_copy:
		p_refined = MergeNodesSubset(graph, p_refined, community_of_main_partition)

	deterior = initial_modularity - modularity(graph, p_refined)
	return sorted(p_refined), deterior


def E_function(graph, comm1, comm2):
	"""This function implements the E(.,.) defined in the paper, the number of links that have one end in comm1 and the other end in comm2.
	
	Arguments:
		graph {[nx.Graph]} -- [the network that will be processed]
		comm1 {[list]} -- [first community that is considered]
		comm2 {[type]} -- [second community that is considered]
	
	Returns:
		[int] -- [the number of links that have one end in C and the other end in D]
	"""
	# to make sure the given communities are in form of lists, even they have only one item.
	if type(comm1) == int:
		comm1 = [comm1]
	if type(comm2) == int:
		comm2 = [comm2]
	
	total_links = 0.0
	for node1, node2 in product(comm1, comm2):
		if graph.has_edge(node1, node2):
			total_links += graph[node1][node2].get('weight', 1.0)

	return total_links


def delta_h(graph, node, old_comm_index, new_comm_index, partition):
	"""computes the improvement of quality function after "node" moves to "community".
	
	Arguments:
		graph {[nx.Graph]} -- [network that will be processed]
		node {[int]} -- [a node in the network]
		community {[list]} -- [a list containing nodes]
		partition {[list]} -- [contains lists that each one represents a community in the network]
	
	Returns:
		[float] -- [the improvement after "node" moves to "community"]
	"""
	if old_comm_index == new_comm_index:
		return 0.0

	increment = 0.0
	decrement = 0.0

	for neigh in graph.neighbors(node):
		if neigh in partition[new_comm_index]:
			increment += graph[node][neigh].get('weight', 1.0)
		elif neigh in partition[old_comm_index]:
			if neigh == node:	# needs to be removed
				continue
			decrement += graph[node][neigh].get('weight', 1.0)

	comm1_size = RecursiveSize(graph, partition[old_comm_index])
	comm2_size = RecursiveSize(graph, partition[new_comm_index])
	node_size = RecursiveSize(graph, [node])

	m = max(1, graph.size(weight='weight'))
	return increment - decrement + (__GAMMA / (2.0 * m)) * node_size * (comm1_size - node_size - comm2_size)


def pick_C_prime_randomly(probabilities):
	"""This function uses the improvement values to choose a community that causes an improvement, but
	not the best improvement necessarily.
	
	Arguments:
		choices {[list]} -- [the list of communities that we want to select randomly one item out of it]
		probabilities {[list]} -- [originally it is the delta_h improvement list, that we consider as probabilities]
	
	Returns:
		[list] -- [the randomly selected community]
	"""
	# if the improvement value is negative, it will not be chosen.
	for i in range(len(probabilities)):
		if probabilities[i] < 0.0:
			probabilities[i] = 0.0
		else:
			try:
				probabilities[i] = exp(probabilities[i] / __THETA)
			except:
				probabilities[i] = __PREDEFINED_INF


	sum_val = sum(probabilities)
	
	# if all movements result in a worse clustering, no new community will be returned.
	if sum_val == 0.0:
		return None

	norm_probabilities = [prob_val / sum_val for prob_val in probabilities]
	return choice(list(range(len(norm_probabilities))), 1, p=norm_probabilities)[0]


def move_node_to_new_comm(partition, node, old_comm_index, new_comm_index):
	"""moves the node to new community "C_prime" by adding the node to it and removing the node
	from the original community. removes the former community if it becomes empty after the move.
	
	Arguments:
		partition {[list]} -- [contains lists that each one represents a community in the network]
		C_prime {[list]} -- [the community that "node" is going to merge in]
		node {[int]} -- [represents a node of the network]
	"""
	if old_comm_index != new_comm_index:
		partition[old_comm_index].remove(node)
		partition[new_comm_index].append(node)

		# if not partition[old_comm_index]:
		# 	partition.pop(old_comm_index)


def MergeNodesSubset(graph, partition, subset):
	"""This function merges the singleton communities in the refined partition. The notations in this function are 
	based on the notation of the paper.
	"""
	R = list()
	remaining_nodes = list(subset)
	for v in subset:
		remaining_nodes.remove(v)
		Kv = graph.degree[v]
		well_connected_threshold = __GAMMA * Kv * (RecursiveSize(graph, subset) - Kv)
		if E_function(graph, v, remaining_nodes) >= well_connected_threshold:
			R.append(v)
		remaining_nodes.append(v)

	T = list()
	quality_improvs = list()
	for v in R:
		old_comm_index = FindNodeInPartition(partition, v)
		if old_comm_index != -1 and len(partition[old_comm_index]) == 1:		# seems iffy, but it really is not :)
			for comm_index in range(len(partition)):
				if set(partition[comm_index]).issubset(set(subset)):
					size_C = RecursiveSize(graph, partition[comm_index])
					well_connected_threshold = __GAMMA * size_C * (RecursiveSize(graph, subset) - size_C)
					if E_function(graph, partition[comm_index], list(set(subset) - set(partition[comm_index]))) >= well_connected_threshold:
						T.append(comm_index)

			for new_comm_index in T:
				quality_improvs.append(delta_h(graph, v, old_comm_index, new_comm_index, partition))
			
			if quality_improvs:
				new_random_index = pick_C_prime_randomly(quality_improvs)
				if new_random_index != None:
					move_node_to_new_comm(partition, v, old_comm_index, T[new_random_index])
					old_comm_index = T[new_random_index]
			
		T.clear()
		quality_improvs.clear()

	i = 0
	while i < len(partition):
		if partition[i] == []:
			partition.pop(i)
			i -= 1
		i += 1

	return partition


def MoveNodesFast(graph, partition):
	"""This function moves nodes between communities until no further move is possible according to local smart move approach.
	The notations in this function are based on the notation of the paper.
	"""
	partition_copy = deepcopy(partition)
	Q = list(graph.nodes())
	shuffle(Q)
	total_improv = 0.0
	improvements = list()
	N = list()
	while True:
		v = Q.pop()
		old_comm_index = FindNodeInPartition(partition_copy, v)
		if partition_copy[-1] != []:
			partition_copy.append([])
		for new_comm_index in range(len(partition_copy)):
			improvements.append(delta_h(graph, v, old_comm_index, new_comm_index, partition_copy))

		max_improvement_index = argmax(improvements)        
		if improvements[max_improvement_index] > 0.0 and old_comm_index != max_improvement_index:
			total_improv +=  improvements[max_improvement_index]
			move_node_to_new_comm(partition_copy, v, old_comm_index, max_improvement_index)
			old_comm_index = max_improvement_index
			N = list(set(graph.neighbors(v)) - set(partition_copy[max_improvement_index]))
			N = list(set(N) - set(Q))

			for node_to_add in N:
				Q.append(node_to_add)
		
		improvements.clear()
		if not Q:
			break
	
	# to remove any empty communities in the partition
	i = 0
	while i < len(partition_copy):
		if partition_copy[i] == []:
			partition_copy.pop(i)
			i -= 1
		i += 1

	return partition_copy, total_improv


def repartition(graph):
	"""nodes of the new graph are the communities of partition_refined. however, the communities 
	of the new graph are the communities of partition. This function determines the initial communities
	of the nodes of the new graph.
	"""
	partition = dict()
	for node in graph.nodes(data=True):
		comm = node[1].get('comm', -1)
		partition[comm] = partition.get(comm, [])
		partition[comm].append(node[0])

	return [partition[key] for key in partition]


def Leiden(graph, partition):
	"""This is the main function that runs the Leidein algorithm on the network, given every node in its own community at the beginning.
	"""
	counter = 0
	current_modularity = 0.0
	current_deterior = 0.0
	current_improv = 0.0

	while True:
		counter += 1

		partition, current_improv = MoveNodesFast(graph, partition)
		current_modularity += (current_improv - current_deterior)
		done = len(partition) == graph.number_of_nodes()

		print('current_modularity at step ', counter, ' = ', current_modularity, end ='\t')
		print('len of p =', len(partition), ' |V(G)| =', graph.number_of_nodes(), 'counter =', counter)

		done = done or (current_improv - current_deterior < __IMPROV_THRESHOLD)
		done = done or (counter > __ITER_THRESHOLD)

		if not done:
			p_refined, current_deterior = RefinePartition(graph, partition)
			print(len(p_refined))
			graph = AggregateGraph(graph, p_refined, partition)
			partition = repartition(graph)

		if done:
			break

	return partition


def extract_info_from_hierarchy(hierarchy, n):
	"""This function finds the true community info based on teh hierarchical structure that is made previously in other functions of the code.
	"""
	if n == 0:
		return hierarchy[n]
	
	lower_level = hierarchy[n-1]
	partition = hierarchy[n]
	
	for comm in partition:
		index_counter = len(comm) - 1
		while index_counter > -1:
			index_in_lower_level = comm[index_counter]
			comm.pop(index_counter)
			to_be_replaced = lower_level[index_in_lower_level]
			comm.extend(to_be_replaced)
			index_counter -= 1

	hierarchy.pop(n - 1)
	return extract_info_from_hierarchy(hierarchy, n-1)


def __main():
	args = utils.create_argument_parser()
	g = utils.load_graph(args.dataset, args.w, report=False)

	initial_partition = SingletonPartition(g)
	L = Leiden(g, initial_partition)
	
	print(L)
	# partition = extract_info_from_hierarchy(L, len(L) - 1)
	# for i, p in enumerate(partition):
	# 	print(i, sorted(p), len(p))

if __name__ == "__main__":
	__main()

