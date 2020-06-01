import networkx as nx
from community_status import Status
import utils
import time
from measures import modularity
from measures import NMI


__PASS_MAX = -1
__MIN = 0.0000001


def best_partition(graph, partition=None, weight='weight', weighted=False, randomize=False):
	"""
	find communities in graph based on SIWO algorithm
	"""


	# ---------------------Step1. Pre-Processing------------------------------
	graph = graph.copy()
	dangles = pre_processing(graph, weighted, weight)

	# ---------------------Step2. Optimizing SIWO-----------------------------

	dendo = generate_dendrogram(graph, inc_SIWO, partition, weight, randomize)
	partition = partition_at_level(dendo, len(dendo) - 1)

	# -----------------Step3. Qualified Community detection-------------

	graph_cpy = graph.copy()
	remove_singles(graph, partition)  # temporarily remove Lone communities

	remove_weights(graph)
	ngraph = induced_graph(partition, graph, weight)

	dendo = generate_dendrogram(ngraph, inc_out_degree, None, weight, randomize)

	dendoN = [partition]
	dendoN.extend(dendo)
	partition = partition_at_level(dendoN, len(dendoN) - 1)

	# ---------------------Step4. Post-Processing-----------------------------

	graph = graph_cpy
	add_singles(partition, graph)
	partition = renumber(partition)

	add_dangles(graph, partition, dangles)

	return partition


def pre_processing(graph, weighted, weight):
	if weighted:
		ccoef = compute_edge_strength(graph, weight='strength')
		normalize_weights(graph, ccoef)
		combine_weights(graph)
	else:
		compute_edge_strength(graph, weight)


	# for n1, n2, edge in graph.edges(data=True):
	# 	print(n1, n2, edge)
	# exit(0)

	all_dangles = []
	dangles = remove_dangles(graph, list(graph.nodes()))
	while(dangles):
		all_dangles.append(dangles)
		dangles = remove_dangles(graph, dangles.values())
	return all_dangles


def compute_edge_strength(graph_base, weight="weight"):
	nodes = list(graph_base.nodes())
	nb_vertices = len(nodes)
	mutuals, max_mutuals, total_mutuals = get_mutuals(graph_base)
	ccoef = clustering_coef(graph_base, total_mutuals)
	for i in range(nb_vertices):
		neighbors = list(graph_base.neighbors(nodes[i]))
		max_mutual = max_mutuals.get(nodes[i])

		for neigh in neighbors:
			if (ccoef.get(nodes[i]) < ccoef.get(neigh)):
				continue

			cur_mutuals = mutuals.get(nodes[i])
			# length of each bin
			binlen = 2/(float)(max_mutual+1)
			# min point of the bin
			minpoint = -1 + (cur_mutuals.get(neigh) * binlen)
			# max point of the bin
			maxpoint = minpoint + binlen
			# average point of the bin
			avgpoint = (minpoint + maxpoint)/2.0
			w = avgpoint

			graph_base.add_edge(nodes[i], neigh, weight= w)

	return ccoef


def get_mutuals(graph_base):
	mutuals = {}
	max_mutuals = {}
	total_mutuals = {}
	nodes = list(graph_base.nodes())
	nb_vertices = len(nodes)

	for i in range(nb_vertices):
		mutuals[nodes[i]] = {}
		max_mutuals[nodes[i]] = -1
		total_mutuals[nodes[i]] = 0.0

	for i in range(nb_vertices):
		neighbors = list(graph_base.neighbors(nodes[i]))
		for neigh in neighbors:
			if neigh in mutuals[nodes[i]]:
				continue
			cur_mutual = shared_neighbors_cnt(graph_base, nodes[i], neigh)
			mutuals[nodes[i]][neigh] = cur_mutual
			mutuals[neigh][nodes[i]] = cur_mutual

			total_mutuals[nodes[i]] = total_mutuals[nodes[i]] + cur_mutual
			total_mutuals[neigh] = total_mutuals[neigh] + cur_mutual

			if cur_mutual > max_mutuals[nodes[i]]:
				max_mutuals[nodes[i]] = cur_mutual
			if cur_mutual > max_mutuals[neigh]:
				max_mutuals[neigh] = cur_mutual

		total_mutuals[nodes[i]] = total_mutuals[nodes[i]]/2.0

	return mutuals, max_mutuals, total_mutuals


def shared_neighbors_cnt(graph_base, u, v):
	shared = 0
	if(graph_base.degree(u) > graph_base.degree(v)):
		tmp = u
		u = v
		v = tmp
	neighbors_u = graph_base[u]
	neighbors_v = graph_base[v]
	for n1 in neighbors_u:
		if n1 in neighbors_v:
			shared = shared+1
	return shared


def clustering_coef(graph_base, total_mutuals):
	clust_coef = {}
	nodes = list(graph_base.nodes())
	for node in nodes:
		total_mutual = total_mutuals[node]
		deg = graph_base.degree(node)
		possible_tri = (deg * (deg-1)) / 2
		if possible_tri == 0:
			clust_coef[node] = 1
		else:
			clust_coef[node] = total_mutual/float(possible_tri)
	return clust_coef


def normalize_weights(graph, ccoef, weight="weight"):
	nodes = list(graph.nodes())
	nb_vertices = len(nodes)
	for i in range(nb_vertices):
		node = nodes[i]
		max_weight = -float('Inf')
		min_weight = float('Inf')
		for neighbor, datas in graph[node].items():  # compute the max weight for each node
			edge_weight = datas.get(weight, 1)
			if edge_weight > max_weight:
				max_weight = edge_weight
			if edge_weight < min_weight:
				min_weight = edge_weight
		for neighbor, datas in graph[node].items():
			if (ccoef.get(nodes[i]) < ccoef.get(neighbor)):
				continue
			# scale edge weight to [0,1]
			edge_weight = datas.get(weight, 1)
			if max_weight == min_weight:
				w = 1.0/len(list(graph.neighbors(node)))
			else:
				w = float(edge_weight - min_weight)/(max_weight-min_weight)

			graph.add_edge(nodes[i], neighbor, w_norm= w)


def combine_weights(graph, weight="weight"):
	nodes = list(graph.nodes())
	nb_vertices = len(nodes)
	for i in range(nb_vertices):
		node = nodes[i]
		for neighbor, datas in graph[node].items():  # combine weight and strength of each edge
			if neighbor < node:
				continue
			edge_weight = datas.get("w_norm", 1)
			edge_strength = datas.get("strength", 0)
			w = ((edge_strength+1.0)/2) + edge_weight-1
			graph.add_edge(nodes[i], neighbor, weight= w)


def remove_dangles(graph, candidates):
	dangles = {}
	for node in candidates:
		deg = graph.degree(node)
		if deg == 1:
			neighbor = list(graph.neighbors(node))[0]
			dangles[node] = neighbor
		else:
			continue
	to_be_rmvd = []
	for node in dangles:
		if graph.degree(node) == 0:
			to_be_rmvd.append(node)
			continue
		graph.remove_node(node)
	for node in to_be_rmvd:
		del dangles[node]
	return dangles


def generate_dendrogram(graph, inc_fnc, part_init=None, weight='weight', randomize=False):
	"""Find communities in the graph and return the associated dendrogram
	A dendrogram is a tree and each level is a partition of the graph nodes.
	Level 0 is the first partition, which contains the smallest communities,
	and the best is len(dendrogram) - 1. The higher the level is, the bigger
	are the communities
	"""

	# special case, when there is no link
	# the best partition is everyone in its community
	if graph.number_of_edges() == 0:
		part = dict([])
		for node in list(graph.nodes()):
			part[node] = node
		return [part]

	current_graph = graph.copy()
	status = Status()
	status.init(current_graph, weight, part_init)
	status_list = list()

	changed = one_level(current_graph, inc_fnc, status, weight, randomize)
	partition = renumber(status.node2com)

	status_list.append(partition)
	current_graph = induced_graph(partition, current_graph, weight)
	status.init(current_graph, weight)

	while False:
		changed = one_level(current_graph, inc_fnc, status, weight, randomize)
		if not changed:
			break
		partition = renumber(status.node2com)
		status_list.append(partition)

		current_graph = induced_graph(partition, current_graph, weight)
		status.init(current_graph, weight)

	return status_list[:]


def one_level(graph, inc_fnc, status, weight_key, randomize=False):
	"""Compute one level of communities
	"""
	modified = True
	changed = False
	while modified:
		print('***')
		modified = False

		for node in randomly(graph.nodes(), randomize):

			com_node = status.node2com[node]
			neigh_communities = neighcom(node, graph, status, weight_key)

			remove_from_community(node, com_node, neigh_communities.get(com_node, 0.), status)
			best_com = com_node
			best_increase = 0.0
			for com, dnc in randomly(neigh_communities.items(), randomize):
				incr = inc_fnc(graph, status, neigh_communities, com, com_node, node)
				if incr > best_increase:
					best_increase = incr
					best_com = com

			insert_to_community(node, best_com, neigh_communities.get(best_com, 0.), status)
			if best_com != com_node:
				modified = True
				changed = True
		
	return changed


def randomly(seq, randomize):
	""" Convert sequence or iterable to an iterable in random order if
	randomize """
	seq = sorted(seq)
	if randomize:
		random.shuffle(seq)
	return seq


def neighcom(node, graph, status, weight_key):
	"""
	Compute the communities in the neighborhood of node in the graph given
	with the decomposition node2com
	"""
	weights = {}
	for neighbor, datas in sorted(graph[node].items()):
		if neighbor != node:
			edge_weight = datas.get(weight_key, 1)
			neighborcom = status.node2com[neighbor]
			weights[neighborcom] = weights.get(neighborcom, 0) + edge_weight

	return weights


def remove_from_community(node, com, weight, status):
	""" Remove node from community com and modify status"""
	status.degrees[com] = (status.degrees.get(com, 0.) - status.gdegrees.get(node, 0.))
	status.internals[com] = float(status.internals.get(com, 0.) - weight - status.loops.get(node, 0.))
	status.node2com[node] = -1
	status.com2nodes[com].remove(node)


def insert_to_community(node, com, weight, status):
	""" Insert node into community and modify status"""
	status.node2com[node] = com
	status.com2nodes.setdefault(com, set()).add(node)
	status.degrees[com] = (status.degrees.get(com, 0.) + status.gdegrees.get(node, 0.))
	status.internals[com] = float(status.internals.get(com, 0.) + weight + status.loops.get(node, 0.))


def renumber(dictionary):
	"""Renumber the values of the dictionary from 0 to n
	"""
	count = 0
	ret = dictionary.copy()
	new_values = dict([])

	for key in dictionary.keys():
		value = dictionary[key]
		new_value = new_values.get(value, -1)
		if new_value == -1:
			new_values[value] = count
			new_value = count
			count += 1
		ret[key] = new_value

	return ret


def induced_graph(partition, graph, weight="weight"):
	"""Produce the graph where nodes are the communities
	there is a link of weight w between communities if the sum of the weights
	of the links between their elements is w
	"""
	ret = nx.Graph()
	ret.add_nodes_from(partition.values())

	for node1, node2, datas in graph.edges(data=True):
		edge_weight = datas.get(weight, 1)
		com1 = partition[node1]
		com2 = partition[node2]
		w_prec = ret.get_edge_data(com1, com2, {weight: 0}).get(weight, 1)
		ret.add_edge(com1, com2, weight= w_prec + edge_weight)

	return ret


def partition_at_level(dendrogram, level):
	"""Return the partition of the nodes at the given level
	A dendrogram is a tree and each level is a partition of the graph nodes.
	Level 0 is the first partition, which contains the smallest communities,
	and the best is len(dendrogram) - 1.
	The higher the level is, the bigger are the communities
	"""
	partition = dendrogram[0].copy()
	for index in range(1, level + 1):
		for node, community in partition.items():
			partition[node] = dendrogram[index][community]
	return partition


def remove_singles(graph, partition):
	singles, com2node = get_single_nodes(partition)
	for node in singles:
		graph.remove_node(node)


def get_single_nodes(partition):
	com2node = {}
	single_nodes = []
	nodes = partition.keys()
	coms = set(partition.values())
	for com in coms:
		com2node[com] = []

	for node in nodes:
		com = partition[node]
		com2node[com].append(node)
	for com in coms:
		if len(com2node[com]) == 1:
			single_nodes.append(com2node[com][0])
	return single_nodes, com2node


def remove_weights(graph, weight="weight"):
	for u, v in graph.edges():
		graph.add_edge(u, v, weight= 1)


def add_singles(partition, graph_base):
	single_nodes, com2node = get_single_nodes(partition)
	single_nodes_cpy = list(single_nodes)

	modified = True
	to_be_rmvd = []

	if len(single_nodes) == 0:
		return partition

	while(modified):
		modified = False

		for node in to_be_rmvd:
			single_nodes.remove(node)

		to_be_rmvd = []
		for node in single_nodes_cpy:
			com_node = partition[node]
			best_com = com_node
			max_com_size = 0
			neighbors = list(graph_base.neighbors(node))
			neigh_com_size = {}  # map between neighbor community and its community size

			for neigh in neighbors:
				if neigh in single_nodes:
					continue

				neighcom = partition[neigh]
				if neighcom not in neigh_com_size:
					neigh_com_size[neighcom] = 1
				else:
					neigh_com_size[neighcom] += 1

			if len(neigh_com_size) == 0:
				continue
			max_com_size = max(neigh_com_size.values())
			min_com_size = min(neigh_com_size.values())
			if(max_com_size == 1):
				# choose neighbor node with largest degree centrality
				max_ds = 0
				for neigh in neighbors:
					if neigh in single_nodes:
						continue
					neighcom = partition[neigh]
					ds = degree_centrality(neigh, len(com2node[neighcom]), partition, graph_base)
					if ds >= max_ds:
						max_ds = ds
						best_com = partition[neigh]
			else:
				for com, size in neigh_com_size.items():
					if size == max_com_size:
						best_com = com
						break

			partition[node] = best_com
			if(best_com != com_node):
				if node in single_nodes:
					to_be_rmvd.append(node)
				modified = True
	if len(single_nodes) != 0:
		# we were not able to add some of the single nodes to the graph (for example when the graph is a circle of nodes)
		com = max(com2node.keys())
		for node in single_nodes:
			partition[node] = com  # put each in its own community
			com += 1
	return partition


def degree_centrality(node, comSize, partition, graph):
	neighbors = list(graph.neighbors(node))
	connected = 0.0
	for neigh in neighbors:
		com_neigh = partition[neigh]

		if com_neigh == partition[node]:
			connected += 1
	return connected/float(comSize-1)


def add_dangles(graph, partition, all_dangles):
	i = len(all_dangles)-1
	while i >= 0:
		dangles = all_dangles[i]
		for u, v in dangles.items():
			partition[u] = partition[v]
		i = i-1


def inc_SIWO(graph, status, neigh_communities, com_neigh, com_node, node):
	return neigh_communities.get(com_neigh, 0.) - neigh_communities.get(com_node, 0.)


def inc_out_degree(graph, status, neigh_communities, com_neigh, com_node, node):

	# in_degrees and out_degrees before and after moving node from com_node to com_neigh

	indeg_node_b = status.loops.get(node, 0.) * 2.0
	outdeg_node_b = status.gdegrees[node] - indeg_node_b
	outdeg_node_a = outdeg_node_b - neigh_communities.get(com_neigh, 0.)

	between_links = neigh_communities.values()
	within_links = indeg_node_b / 2.0
	max_between_links = max(between_links)
	if indeg_node_b > outdeg_node_b and within_links > max_between_links:
		return 0
	incr = outdeg_node_b - outdeg_node_a
	return incr




def main():
	start_time = time.time()

	args = utils.create_argument_parser()
	graph = utils.load_graph(args.dataset, args.w)
	partition = best_partition(graph)

	communities = utils.extract_communities(partition)
	# utils.print_comm_info_to_display(communities)
	utils.write_comm_info_to_file(args.output, partition)

	# print('modularity_value =', modularity(graph, communities))
	# print('NMI =', NMI(args.output, partition))

	finish_time = time.time()
	print('\nDone in %.4f seconds.' %(finish_time - start_time))



if __name__ == "__main__":
	main()