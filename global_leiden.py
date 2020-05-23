import argparse
import networkx as nx
import random
from copy import deepcopy
from itertools import product
import numpy as np
import time


__GAMMA = 1.0


def create_argument_parser():
    # This function reads the command line argument and parse them  
	parser = argparse.ArgumentParser(description='The Leiden algorithm', prog="Leiden")
	parser.add_argument('dataset', help='Path to the network files')
	parser.add_argument('output', help='output name')
	parser.add_argument('gamma', help='resolution parameter')
	parser.add_argument('-w', help='weighted graphs', action="store_true")
	parser.add_argument('-r', help='randomize order of nodes', action="store_true")
	return parser.parse_args()


def load_graph(fileaname):
	# load the graph from a file containing list of edges (unweighted)
	graph = nx.Graph()
	with open(fileaname, 'r') as file:
		lines = file.readlines()
		for line in lines:
			nodes = line.split()
			nodes = list(map(int, nodes))
			graph.add_node(nodes[0], size= 1, inner_nodes = [nodes[0]])
			graph.add_node(nodes[1], size= 1, inner_nodes = [nodes[1]])
			graph.add_edge(nodes[0], nodes[1], weight = 1.0)

	return graph


def make_singleton_communities(graph):
	# put each node in a community of size one
	partition = list()
	for node in sorted(graph.nodes()):
		partition.append([node])
	return partition


def find_node_in_partition(partition, node):
	# iterate through communities in the partition to find the right community index
	for com_index in range(len(partition)):
		for node_in_com in partition[com_index]:
			if node == node_in_com:
				return com_index
	# this should not happen, if node is not inside any community, -1 is returned
	print('ERROR IN **find_node_in_partition**')
	return -1


def aggregate_graph(graph, partition):
	# create a new graph based on the given partition in which every node represents community in the given graph
	new_agg = nx.Graph()
	# add nodes to the new aggregated graph
	for com_index in range(len(partition)):
		node_size = 0
		inner_nodes = list()
		for node in partition[com_index]:
			node_size += get_size_node(graph, node)
			inner_nodes.extend(graph.nodes[node].get('inner_nodes', []))
		new_agg.add_node(com_index, size = node_size, inner_nodes = inner_nodes)

	# add edges to the new aggregated graph
	for n1, n2 , edge in graph.edges(data=True):
		com_index_node1 = find_node_in_partition(partition, n1)
		com_index_node2 = find_node_in_partition(partition, n2)
		w = edge.get('weight', 1.0)
		if new_agg.has_edge(com_index_node1, com_index_node2):
			w += new_agg[com_index_node1][com_index_node2].get('weight', 1.0)
		new_agg.add_edge(com_index_node1, com_index_node2, weight = w)
	return new_agg


def get_size_node(graph, node):
	# get the number of actual node in a given relative node
	node_info = graph.nodes[node]
	inner_nodes = node_info.get('inner_nodes', [])
	return max(1, len(inner_nodes))


def get_size_community(graph, community):
	# get the number of actual nodes in the given community
	com_size = 0
	for node in community:
		com_size += get_size_node(graph, node)
	return com_size


def get_size_partition(graph, partition):
	# get the number of actual nodes in the given partition
	partition_size = 0
	for community in partition:
		partition_size += get_size_community(graph, community)
	return partition_size


# def modularity(graph, partition):
# 	# computes the modularity for the graph considering the given partition
# 	sum_weights = 0.0
# 	for n1, n2, edge in graph.edges(data=True):
# 		com_index1 = find_node_in_partition(partition, n1)
# 		com_index2 = find_node_in_partition(partition, n2)
# 		if com_index1 == com_index2:
# 			sum_weights += edge.get('weight', 1.0)

# 	sum_com_sizes = 0.0
# 	for com_index in range(len(partition)):
# 		com_size = get_size_community(graph, partition[com_index])
# 		sum_com_sizes += com_size * (com_size - 1) / 2.0

# 	m = graph.size(weight='weight')
# 	return (sum_weights - __GAMMA * sum_com_sizes) / (2.0 * m)


def delta_modularity(graph, node, partition, curr_com_index, next_com_index):
	# compute how much modularity changes if a node moves from the current community to the new community
	mod = 0.0

	if curr_com_index == next_com_index:
		return mod

	for neigh in graph.neighbors(node):
		if neigh == node:
			continue
		if neigh in partition[curr_com_index]:
			mod -= graph[node][neigh].get('weight', 1.0)
		elif neigh in partition[next_com_index]:
			mod += graph[node][neigh].get('weight', 1.0)
	
	node_size = get_size_node(graph, node)
	curr_com_size = get_size_community(graph, partition[curr_com_index])
	next_com_size = get_size_community(graph, partition[next_com_index])
	m = graph.size(weight='weight')

	mod += __GAMMA * node_size * (curr_com_size - node_size - next_com_size)
	return mod / (2.0 * m)


def move_node(node, partition, old_com_index, new_com_index):
	# move the given node inside the partition from the old community to the new community
	if old_com_index == new_com_index:
		print('NO MOVE IS DONE')
		return new_com_index

	partition[new_com_index].append(node)
	partition[old_com_index].remove(node)
	return new_com_index


def find_neighbor_coms(graph, node, partition):
	# find the indices of the communities that are connected to the given node (current community is excluded)
	curr_com_index = find_node_in_partition(partition, node)
	neighbor_com_indices = list()
	for neigh in graph.neighbors(node):
		com_index_of_neigh = find_node_in_partition(partition, neigh)
		neighbor_com_indices.append(com_index_of_neigh)
	neighbor_com_indices = list(set(neighbor_com_indices) - set([curr_com_index]))
	return neighbor_com_indices


# def find_random_for_tiebreak(partition, indices):
# 	max_index = indices[0]
# 	for i in range(1, len(indices)):
# 		if len(partition[indices[i]]) > len(partition[max_index]):
# 			max_index = indices[i]
# 	return max_index


def move_nodes_fast(graph, partition, report=False, mode='MOVE'):
	# runs the first step of Leiden: moving nodes to neighboring communities
	queue_of_nodes = list()
	if mode == 'MOVE':
		queue_of_nodes = [(x, 'MOVE') for x in list(graph.nodes())]
	elif mode == 'MERGE':
		for com_index in range(len(partition)):
			if len(partition[com_index]) == 1:
				queue_of_nodes.append((partition[com_index][0], 'MERGE'))

	# random.shuffle(queue_of_nodes)
	
	# add one empty community for when moving to an empty community is required
	partition.append([])

	while queue_of_nodes:
		# determine the node to be investigated
		(node, state) = queue_of_nodes.pop(0)
		curr_com_index = find_node_in_partition(partition, node)

		delta_mod = dict()
		# evaluate how much the quality of the partition changes, if the node moves to a neighboring community
		neigh_coms = find_neighbor_coms(graph, node, partition)
		for com_index in neigh_coms:
			delta_mod[com_index] = delta_modularity(graph, node, partition, curr_com_index, com_index)

		# if remains in its current community, means no movement, and there is no improvement for it
		delta_mod[curr_com_index] = 0.0

		# moving to an empty community should be evaluated too if the current community has more than one node
		# since |C| == |V|, there always is at least one empty community if the below condition is satisfied.
		if len(partition[curr_com_index]) != 1 and mode == 'MOVE':
			try:
				empty_com_index = partition.index([])
			except Exception as e:
				partition.append([])
				empty_com_index = len(partition) - 1
			delta_mod[empty_com_index] = delta_modularity(graph, node, partition, curr_com_index, empty_com_index)

		if delta_mod:
			max_delta_mod_index = curr_com_index
			for index, value in delta_mod.items():
				if value > delta_mod[max_delta_mod_index]:
					max_delta_mod_index = index

			if curr_com_index != max_delta_mod_index:
				curr_com_index = move_node(node, partition, curr_com_index, max_delta_mod_index)

				# adding old neighbors of the node to the queue if they are not there, to be investigated again
				old_neighbors = [x for x in list(graph.neighbors(node)) if not x in partition[max_delta_mod_index]]
				for neigh in [x for x in old_neighbors if not (x, 'MOVE') in queue_of_nodes]:
					queue_of_nodes.append((neigh, 'MOVE'))

	remove_empty_coms(partition)
	return partition


def merge_nodes_fast(graph, partition, report=False):
	# runs the first step of Leiden: merging nodes in communities of size one to neighboring communities
	return move_nodes_fast(graph, partition, report=report, mode='MERGE')


def remove_empty_coms(partition):
	# remove any empty element in the partition
	i = 0
	while i < len(partition):
		if len(partition[i]) == 0:
			partition.pop(i)
			i -= 1
		i += 1
	return partition 


def move_nodes_refine(graph, partition, community, report=False, mode='MOVE'):
	# R is the set of nodes that are well connected inside the given community
	queue_of_nodes = list()
	if mode == 'MOVE':
		queue_of_nodes = [(x, 'MOVE') for x in list(graph.nodes())]
	elif mode == 'MERGE':
		for com_index in range(len(partition)):
			if len(partition[com_index]) == 1:
				queue_of_nodes.append((partition[com_index][0], 'MERGE'))

	# random.shuffle(queue_of_nodes)
	
	# add one empty community for when moving to an empty community is required
	partition.append([])

	while queue_of_nodes:
		# determine the node to be investigated
		(node, state) = queue_of_nodes.pop(0)
		curr_com_index = find_node_in_partition(partition, node)

		delta_mod = dict()
		# evaluate how much the quality of the partition changes, if the node moves to a neighboring community
		neigh_coms = find_neighbor_coms(graph, node, partition)
		for com_index in neigh_coms:
			delta_mod[com_index] = delta_modularity(graph, node, partition, curr_com_index, com_index)

		# if remains in its current community, means no movement, and there is no improvement for it
		delta_mod[curr_com_index] = 0.0

		# moving to an empty community should be evaluated too if the current community has more than one node
		# since |C| == |V|, there always is at least one empty community if the below condition is satisfied.
		if len(partition[curr_com_index]) != 1 and mode == 'MOVE':
			try:
				empty_com_index = partition.index([])
			except Exception as e:
				partition.append([])
				empty_com_index = len(partition) - 1
			delta_mod[empty_com_index] = delta_modularity(graph, node, partition, curr_com_index, empty_com_index)

		if delta_mod:
			max_delta_mod_index = curr_com_index
			for index, value in delta_mod.items():
				if value > delta_mod[max_delta_mod_index]:
					max_delta_mod_index = index

			if curr_com_index != max_delta_mod_index:
				curr_com_index = move_node(node, partition, curr_com_index, max_delta_mod_index)

				# adding old neighbors of the node to the queue if they are not there, to be investigated again
				old_neighbors = [x for x in list(graph.neighbors(node)) if not x in partition[max_delta_mod_index]]
				for neigh in [x for x in old_neighbors if not (x, 'MOVE') in queue_of_nodes]:
					queue_of_nodes.append((neigh, 'MOVE'))

	remove_empty_coms(partition)
	return partition


def merge_nodes_refine(graph, partition, community, report):
	# runs the first step of Leiden: merging nodes in communities of size one to neighboring communities
	return move_nodes_refine(graph, partition, community, report=report, mode='MERGE')


def refine_partition(graph, partition, report=False):
	# refine the original partition to potentially more and smaller communities
	partition_copy = deepcopy(partition)
	graph_copy = deepcopy(graph)
	refined = list()
	for community in partition_copy:
		smaller_graph = graph_copy.subgraph(community)
		splitted = make_singleton_communities(smaller_graph)
		splitted = move_nodes_refine(smaller_graph, splitted, community, report=report)
		splitted = merge_nodes_refine(smaller_graph, splitted, community, report=report)
		refined.extend(splitted)
	return refined


def repartition(graph, refined, partition):
	# determine the initial partition for the given graph based on the original previous partition (before refinement step)
	new_partition = dict()
	for node in graph.nodes():
		node_in_last_partition = refined[node][0]
		com_index = find_node_in_partition(partition, node_in_last_partition)
		new_partition[com_index] = new_partition.get(com_index, [])
		new_partition[com_index].append(node)
	return [new_partition[k] for k in new_partition]


def run_leiden(graph, partition, report=False, num_iter_limit=2):
	# run the Leiden algorithm on the graph based on the given initial partition
	num_iter = 0
	len_partition = -1

	print('Running Leiden alg using resolution_parameter =', __GAMMA)

	while True:
		num_iter += 1
		if report:
			print('num_iter =', num_iter)

		partition = move_nodes_fast(graph, partition, report=report)
		partition = merge_nodes_fast(graph, partition, report=report)

		if report:
			print('partition =', partition, len(partition))

		if len_partition == graph.number_of_nodes() or num_iter == num_iter_limit:
			break

		refined = refine_partition(graph, partition, report=report)

		if report:
			print('refined =', refined, len(refined))

		graph = aggregate_graph(graph, refined)
		partition = repartition(graph, refined, partition)
		len_partition = len(partition)

	print('number of iteration(s) =', num_iter)
	return partition, graph


def extract_final_communities_out_of_partition(partition, graph):
	# extract the true communities out of the partition
	communities = {i:[] for i in range(len(partition))}
	for i in communities:
		for node in partition[i]:
			communities[i].extend(graph.nodes[node].get('inner_nodes', []))
		communities[i] = sorted(communities[i])
	sorted_indices = sorted(communities, key=lambda k:len(communities[k]), reverse=True)
	return {index:communities[index] for index in sorted_indices}


def set_resolution_parameter(gamma):
	# set the resolution parameter which is a global variable in this code named "__GAMMA"
	global __GAMMA
	__GAMMA = gamma


def __main():
	args = create_argument_parser()

	start_time = time.time()
	set_resolution_parameter(float(args.gamma))
	graph = load_graph(args.dataset)

	partition = make_singleton_communities(graph)
	partition, graph = run_leiden(graph, partition, report=False)
	communities = extract_final_communities_out_of_partition(partition, graph)
	
	for com_index, nodes in communities.items():
		print(com_index, ':', nodes, '\t', len(nodes))

	finish_time = time.time()
	print('\nDone in %.4f seconds.' %(finish_time - start_time))


if __name__ == "__main__":
	__main()

