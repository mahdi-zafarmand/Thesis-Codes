import networkx as nx
import utils
import time
from math import sqrt


def compute_m(graph, community, new_node):
	"""computes the local modularity M if the new node joins the community

	Args:
		graph ([nx.Graph]): [the given network]
		community ([list]): [list of nodes that are already discovered]
		new_node ([int]): [the expectied node to join the community]

	Returns:
		[float]: [the modularity R if the new node joins the community]
	"""
	new_community = community + [new_node]
	internal_edges = 0
	external_edges = 0

	for node in new_community:
		for neigh in graph.neighbors(node):
			if neigh in new_community:
				internal_edges += 1
			else:
				external_edges += 1

	return internal_edges / external_edges	


def find_the_best_neighbor(mod_m):
	"""finds the best node from the neighborhood to join the community based on highest modularity M

	Args:
		mod_m ([dict]): [keys: nodes in neighborhood, values: the corresponding modularity M]

	Returns:
		[int]: [a node from the neighborhood]
	"""
	neighbors = list(mod_m.keys())
	best_neigh = neighbors[0]
	for i in range(1, len(neighbors)):
		if mod_m[neighbors[i]] > mod_m[best_neigh]:
			best_neigh = neighbors[i]

	return best_neigh


def update_neighbors(graph, community, the_neighbors, best_neigh, all_communities=[]):
	"""updates the list of neighbors after a best candidate merges into the community.

	Args:
		graph ([nx.Graph]): [the given network]
		community ([list]): [list of already discovered nodes]
		the_neighbors ([list]): [the previous list of neighbors]
		best_neigh ([int]): [the former best candidate that joined the community]
		all_communities ([list]): [a list of all detected communities]

	Returns:
		[list]: [the updated list of neighbors after the best candidate joins the community]
	"""
	neighbors_of_best_neigh = list(set(graph.neighbors(best_neigh)) - set(community))
	if neighbors_of_best_neigh != []:
		the_neighbors.extend(neighbors_of_best_neigh)
		the_neighbors = list(set(the_neighbors))
	the_neighbors.remove(best_neigh)

	if len(all_communities) > 0:
		for node in [dummy_node for community in all_communities for dummy_node in community]:
			if node in the_neighbors:
				the_neighbors.remove(node)

	return the_neighbors


def update_boundaries(graph, community, boundary, best_neigh):
	""" adds the best_neigh to the boundary if necessary (if the best_neigh is only
	connected to the core, it won't be added to boundaries)

	Args:
		graph ([nx.Graph]): [the given network]
		community ([list]): [list of already discovered nodes]
		boundary ([list]): [group of nodes in the community that has at least one neighbor outside of it]
		best_neigh ([int]): [the former best candidate that joined the community]

	Returns:
		[list]: [the updated boundary list]
	"""
	not_add_to_boundary = True
	for neigh in graph.neighbors(best_neigh):
		if not neigh in community:
			not_add_to_boundary = False
			break

	if not not_add_to_boundary:
		boundary.append(best_neigh)

	# find any node in boundary that should be updated to go to core
	nodes_to_remove = list()
	for node in boundary:
		should_update = True
		for neigh in graph.neighbors(node):
			if not neigh in community:
				should_update = False
				break
		if should_update:
			nodes_to_remove.append(node)

	# remove any node in boundary that should be updated to go to core
	for node in nodes_to_remove:
		boundary.remove(node)

	return boundary


def community_search(graph, intended_node):
	"""	searches for all nodes in the network that should be in the same community as the $intended_node.

	Args:
		graph ([nx.Graph]): [the given network]
		intended_node ([int]): [a node in the network that starts the community expansion]

	Returns:
		[list]: [the discovered community]
	"""
	community = list()
	boundary = list()
	neighbors = list()

	M_measure = 0.0
	community.append(intended_node)
	boundary.append(intended_node)
	neighbors = list(graph.neighbors(intended_node))

	while len(community) < graph.number_of_nodes():
		# compute the improvement caused by any node in the neighborhood
		mod_m = dict()
		for neigh in neighbors:
			mod_m[neigh] = compute_m(graph, community, neigh)

		# find the neighbor with the highest improvement to the modularity R
		best_neigh = find_the_best_neighbor(mod_m)

		# continue expanding the community only if modularity R increases
		if mod_m[best_neigh] < M_measure:
			break

		# update modularity R and add the best neighbor to the community
		M_measure = mod_m[best_neigh]
		community.append(best_neigh)

		# update the neighborhood list
		neighbors = update_neighbors(graph, community, neighbors, best_neigh)

		# update the boundary list
		boundary = update_boundaries(graph, community, boundary, best_neigh)

	return sorted(community)


def community_detection(graph, ground_truth_file_address):
	"""detects all communities in the network one after the other using only local information and updation.

	Args:
		graph ([nx.Graph]): [the given network]
		ground_truth_file_address ([str]): [filename of the ground-truth information of communities]

	Returns:
		[list]: [all discovered communities]
	"""
	start_time = time.time()
	all_communities = list()
	nodes_discovered = 0

	while nodes_discovered < graph.number_of_nodes():
		intended_node = utils.find_node_highest_degree(graph, all_communities)
		community = list()
		boundary = list()
		neighbors = list()
		M_measure = 0.0

		community.append(intended_node)
		nodes_discovered += 1
		boundary.append(intended_node)
		neighbors = list(graph.neighbors(intended_node))
		for node in [dummy_node for community in all_communities for dummy_node in community]:
			if node in neighbors:
				neighbors.remove(node)

		while len(community) < graph.number_of_nodes():

			# compute the improvement caused by any node in the neighborhood
			mod_m = dict()
			for neigh in neighbors:
				mod_m[neigh] = compute_m(graph, community, neigh)

			if len(mod_m) == 0:
				break

			# find the neighbor with the highest improvement to the modularity R
			best_neigh = find_the_best_neighbor(mod_m)

			# continue expanding the community only if modularity R increases
			if mod_m[best_neigh] < M_measure:
				break

			# update modularity R and add the best neighbor to the community
			M_measure = mod_m[best_neigh]
			community.append(best_neigh)
			nodes_discovered += 1

			# update the neighborhood list
			neighbors = update_neighbors(graph, community, neighbors, best_neigh, all_communities)

			# update the boundary list
			boundary = update_boundaries(graph, community, boundary, best_neigh)

		all_communities.append(community)

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
