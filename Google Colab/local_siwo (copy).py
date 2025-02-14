import networkx as nx
import argparse
import utils
import measures
import time
from math import sqrt
from copy import deepcopy


__MIN = 0.0000001


def update_mutual_dicts(graph, node, mutuals, max_mutuals):
	"""updates the mutual information for every two connected nodes in the network gradually.
	the main reason is after each call, we have full info of $node and a little info of its neighbors (will update later).

	Args:
		graph ([nx.Graph]): [the given network]
		node ([int]): [the main node of the grpah which its info is going to be updated]
		mutuals ([dict]): [keys: nodes of network, values: dicts(keys: neighbors of node, values: number of common neighbor)]
		max_mutuals ([dict]): [keys: nodes of the network, values: highest number of common neighbors between key and its neighbors]
	"""
	neighbors = list(graph.neighbors(node))
	# initilization for nodes that haven't been exposed to this function before
	if not node in mutuals:
		mutuals[node] = {}
		max_mutuals[node] = -1

	# each time, updates info for two nodes, node and one of its neighbors (the neighbor will be updated later too.)
	for neigh in neighbors:
		if neigh in mutuals[node]:
			continue
		if not neigh in mutuals:
			mutuals[neigh] = {}
			max_mutuals[neigh] = -1

		cur_mutual = sum(1 for _ in nx.common_neighbors(graph, node, neigh))
		mutuals[node][neigh] = cur_mutual
		mutuals[neigh][node] = cur_mutual

		if cur_mutual > max_mutuals[node]:
			max_mutuals[node] = cur_mutual		# it will not change forever
		if cur_mutual > max_mutuals[neigh]:
			max_mutuals[neigh] = cur_mutual		# it will be updated when this function is called with (node=neigh)


def assign_local_strength(graph, node, mutuals, max_mutuals):
	"""assigns or updates local strength of the edges incident to $node.

	Args:
		graph ([nx.Graph]): [int]
		node ([int]): [the main node of the grpah which its edges are going to be updated]
		mutuals ([dict]): [keys: nodes of network, values: dicts(keys: neighbors of node, values: number of common neighbor)]
		max_mutuals ([dict]): [keys: nodes of the network, values: highest number of common neighbors between key and its neighbors]
	"""
	update_mutual_dicts(graph, node, mutuals, max_mutuals)
	max_mutual_node = max_mutuals.get(node)

	for neigh in graph.neighbors(node):
		max_mutual_neigh = max_mutuals.get(neigh)
		w = mutuals.get(node).get(neigh)
		try:
			w1 = w / max_mutual_node
		except:
			w1 = 0.0
		try:
			w2 = w / max_mutual_neigh
		except:
			w2 = 0.0
		w = w1 + w2 - 1.0
		graph.add_edge(node, neigh, strength=w)


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


def find_best_new_candidate(graph, community, candidates):
	"""finds the best candidate among all candidates and return it and the improvement is causes.

	Args:
		graph ([nx.Graph]): [the given network]
		community ([list]): [contains nodes that are already discovered to be in the community]
		candidates ([list]): [all direct neighbors of nodes in the community which are not in the community]

	Returns:
		[tuple]: [best candidate and the improvement is causes]
	"""
	sum_strength = dict()
	for candidate in candidates:
		sum_strength[candidate] = 0.0
		for neigh in graph.neighbors(candidate):
			if neigh in community:
				sum_strength[candidate] += graph[candidate][neigh].get('strength', 0.0)

	best_candidate = candidates[0]
	for candidate in candidates:
		if sum_strength[candidate] > sum_strength[best_candidate]:
			best_candidate = candidate

	return best_candidate, sum_strength[best_candidate]


def amend_by_dangles(graph, community, quality):
	"""amends the discovered community by adding every dangling nodes that are connected to one member of the community.

	Args:
		graph ([nx.Graph]): [the given network]
		community ([list]): [contains nodes that are already discovered to be in the community]

	Returns:
		[list]: [all nodes that should be placed in the community]
	"""
	dangles = list()
	for node in community:
		for neigh in graph.neighbors(node):
			if graph.degree[neigh] == 1:
				dangles.append(neigh)
				quality += graph[node][neigh].get('strength', 0.0)

	return sorted(community + dangles), quality


def community_search(graph, intended_node):
	"""searches for all nodes in the network that should be in the same community as the $intended_node.

	Args:
		graph ([nx.Graph]): [the given network]
		intended_node ([int]): [a node in the network that starts the community expansion]

	Returns:
		[tuple]: [the discovered community and its quality based on SIWO Function]
	"""
	mutuals = dict()
	max_mutuals = dict()

	community = list()
	cur_quality = 0.0

	community.append(intended_node)
	assign_local_strength(graph, intended_node, mutuals, max_mutuals)

	while len(community) < graph.number_of_nodes():
		new_candidates = utils.find_all_new_candidates(graph, community)

		if new_candidates == []:	# end the discovery if there is no new candidate to explore
			break

		for candidate in new_candidates:
			assign_local_strength(graph, candidate, mutuals, max_mutuals)

		new_node, delta_quality = find_best_new_candidate(graph, community, new_candidates)

		if delta_quality < __MIN:	# end the discovery if the best candidate cannot improve the quality
			break

		community.append(new_node)
		cur_quality += delta_quality

	community, cur_quality = amend_by_dangles(graph, community, cur_quality)
	return community, cur_quality


def report_performance(graph, communities, ground_truth_file_address):
	"""prints some performace measure of the experiment, inclusing Qmodularity and NMI score of the partition.

	Args:
		graph ([nx.Graph]): [the given network]
		communities ([list]): [list of all discovered communities]
		ground_truth_file_address ([str]): [filename of the ground-truth information of communities]
	"""
	com2nodes = dict()
	for i in range(len(communities)):
		com2nodes[i] = communities[i]

	node2com = dict()
	for com_index, nodes in com2nodes.items():
		for node in nodes:
			node2com[node] = com_index

	print('Modularity = %.3f' %(measures.modularity(graph, com2nodes)), end='\t')
	print('NMI = %.3f' %(measures.NMI(ground_truth_file_address, node2com)))


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


def find_best2(graph, new_node, candidates, candidates_improvs):
	"""finds the best candidate to join the community based on updated improvements

	Args:
		graph ([nx.Graph]): [the given network]
		new_node ([int]): [the most recent node in the community]
		candidates ([list]): [all candidates to be considered]
		candidates_improvs ([dict]): [candidates' improvements]

	Returns:
		[tuple]: [best candidate and its corresponding improvement]
	"""
	for node in candidates:
		if not node in candidates_improvs:
			candidates_improvs[node] = graph[node][new_node].get('strength', 0.0)
		elif graph.has_edge(node, new_node):
			candidates_improvs[node] += graph[node][new_node].get('strength', 0.0)
	
	if new_node in candidates_improvs:
		del candidates_improvs[new_node]
	
	best_candidate = candidates[0]
	for candidate in candidates:
		if candidates_improvs[candidate] > candidates_improvs[best_candidate]:
			best_candidate = candidate

	return best_candidate, candidates_improvs[best_candidate]


def community_detection(graph, ground_truth_file_address):
	"""detects all communities in the network one after the other using only local information and updation.

	Args:
		graph ([nx.Graph]): [the given network]
		ground_truth_file_address ([str]): [filename of the ground-truth information of communities]

	Returns:
		[list]: [all discovered communities]
	"""
	start_time = time.time()
	graph_copy = deepcopy(graph)
	graph_copy = remove_self_loops(graph_copy)

	mutuals = dict()
	max_mutuals = dict()

	nodes_discovered = 0
	all_communities = list()
	fully_analyzed_nodes = set()

	while nodes_discovered < graph.number_of_nodes():
		community = list()
		candidates = list()
		candidates_improvs = dict()

		initial_node = utils.find_node_highest_degree(graph_copy, all_communities)			# among undiscovered nodes
		community.append(initial_node)														# update community
		nodes_discovered += 1																# update number of discovered nodes
		candidates = list(set(graph_copy.neighbors(initial_node)) - set(community))			# update list of candidates

		assign_local_strength(graph_copy, initial_node, mutuals, max_mutuals)				# of the edges that are connected to the initial node
		fully_analyzed_nodes.add(initial_node)												# update list of fully analyzed nodes
		
		while candidates != []:			

			for candidate_node in candidates:												# to the edges that connects the community to its neighborhood
				if not candidate_node in fully_analyzed_nodes:
					assign_local_strength(graph_copy, candidate_node, mutuals, max_mutuals)
					fully_analyzed_nodes.add(candidate_node)

			new_node, delta_quality = find_best2(graph_copy, community[-1], candidates, candidates_improvs)	# find the best candidate to merge into the community

			if delta_quality < __MIN:										# end the discovery if the best candidate cannot improve the quality
				break				

			community.append(new_node)										# adds the best candidates among all neighbors to the community
			nodes_discovered += 1
			candidates.extend(graph_copy.neighbors(new_node))
			candidates = list(set(candidates) - set(community))

		graph_copy.remove_nodes_from(community)
		all_communities.append(community)										# sort is for better representation only!

	all_communities = amend_partition(graph, all_communities)	
	report_performance(graph, all_communities, ground_truth_file_address)

	print('Community Detection Task is done in %.4f seconds.' %(time.time() - start_time))
	return all_communities


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


def community_search_for_all_nodes(graph, ground_truth_file_address):
	"""do community search for all nodes separately, calculate accuracy measures then report the avg. and sd.

	Args:
		graph ([nx.Graph]): [the given network]
		ground_truth_file_address ([str]): [filename of the ground-truth information of communities]
	"""
	start_time = time.time()
	performance_info = dict()
	ground_truth_com2nodes = utils.read_ground_truth(ground_truth_file_address)

	for e, node in enumerate(graph.nodes()):
		performance_info[node] = {'degree': graph.degree[node], 'precision': 0.0, 'recall': 0.0, 'f1-score': 0.0}
		community = community_search(graph, node)[0]
		update_performance_info(node, performance_info, community, ground_truth_com2nodes)

	precision = sum(performance_info[x]['precision'] for x in list(graph.nodes())) / graph.number_of_nodes()
	recall = sum(performance_info[x]['recall'] for x in list(graph.nodes())) / graph.number_of_nodes()
	f1_score = sum(performance_info[x]['f1_score'] for x in list(graph.nodes())) / graph.number_of_nodes()
	sd = sqrt(sum((performance_info[x]['f1_score'] - f1_score) ** 2 for x in list(graph.nodes())) / graph.number_of_nodes())

	print('precision =', precision, '   recall =', recall, '   f1-score =', f1_score, '   sd(fscore) =', sd)
	finish_time = time.time()
	print('Community Search for all nodes is done in %.4f seconds.' %(finish_time - start_time))

