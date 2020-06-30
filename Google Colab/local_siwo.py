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

		# cur_mutual : number of triangles that edge(node, neigh) belongs to
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


def find_best_new_candidate(graph, new_node, candidates, candidates_improvs):
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


def amend_for_community_search(graph, community):
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

	return sorted(set(community + dangles))


def community_search(graph, intended_node):
	"""searches for all nodes in the network that should be in the same community as the $intended_node.

	Args:
		graph ([nx.Graph]): [the given network]
		intended_node ([int]): [a node in the network that starts the community expansion]

	Returns:
		[list]: [the discovered community]
	"""
	start_time = time.time()
	graph_copy = deepcopy(graph)
	graph_copy = utils.remove_self_loops(graph_copy)

	mutuals = dict()
	max_mutuals = dict()
	fully_analyzed_nodes = set()

	community = list()
	community.append(intended_node)
	candidates = list(set(graph_copy.neighbors(intended_node)) - set(community))
	candidates_improvs = dict()

	assign_local_strength(graph_copy, intended_node, mutuals, max_mutuals)
	fully_analyzed_nodes.add(intended_node)

	while len(community) < graph_copy.number_of_nodes() and candidates != []:

		for candidate_node in candidates:
			if not candidate_node in fully_analyzed_nodes:
				assign_local_strength(graph_copy, candidate_node, mutuals, max_mutuals)
				fully_analyzed_nodes.add(candidate_node)

		new_node, delta_quality = find_best_new_candidate(graph_copy, community[-1], candidates, candidates_improvs)

		if delta_quality < __MIN:	# end the discovery if the best candidate cannot improve the quality
			break

		community.append(new_node)
		candidates.extend(graph_copy.neighbors(new_node))
		candidates = list(set(candidates) - set(community))

	# community = amend_for_community_search(graph_copy, community)
	return community


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
	graph_copy = utils.remove_self_loops(graph_copy)

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

			new_node, delta_quality = find_best_new_candidate(graph_copy, community[-1], candidates, candidates_improvs)

			if delta_quality < __MIN:														# end the discovery if the best candidate cannot improve the quality
				break				

			community.append(new_node)														# adds the best candidates among all neighbors to the community
			nodes_discovered += 1
			candidates.extend(graph_copy.neighbors(new_node))
			candidates = list(set(candidates) - set(community))

		graph_copy.remove_nodes_from(community)
		all_communities.append(community)													# sort is for better representation only!

	all_communities = utils.amend_partition(graph, all_communities)	
	report = utils.report_performance(graph, all_communities, ground_truth_file_address)
	print(report, '\tCommunity Detection Task is done in %.4f seconds.' %(time.time() - start_time))
	return all_communities


def community_search_for_all_nodes(graph, ground_truth_file_address):
	"""do community search for all nodes separately, calculate accuracy measures then report the Avg. and SD.

	Args:
		graph ([nx.Graph]): [the given network]
		ground_truth_file_address ([str]): [filename of the ground-truth information of communities]
	"""
	start_time = time.time()
	performance_info = dict()
	ground_truth_com2nodes = utils.read_ground_truth(ground_truth_file_address)

	for e, node in enumerate(graph.nodes()):
		performance_info[node] = {'degree': graph.degree[node], 'precision': 0.0, 'recall': 0.0, 'f1-score': 0.0}
		community = community_search(graph, node)
		utils.update_performance_info(node, performance_info, community, ground_truth_com2nodes)

	precision = sum(performance_info[x]['precision'] for x in list(graph.nodes())) / graph.number_of_nodes()
	recall = sum(performance_info[x]['recall'] for x in list(graph.nodes())) / graph.number_of_nodes()
	f1_score = sum(performance_info[x]['f1_score'] for x in list(graph.nodes())) / graph.number_of_nodes()
	sd = sqrt(sum((performance_info[x]['f1_score'] - f1_score) ** 2 for x in list(graph.nodes())) / graph.number_of_nodes())

	print('precision = %.4f' %precision, end='\t')
	print('recall = %.4f' %recall, end='\t')
	print('f1-score = %.4f' %f1_score, end='\t')
	print('sd(fscore) = %.4f' %sd, end='\t')
	finish_time = time.time()
	print('time = %.4f' %(finish_time - start_time))
