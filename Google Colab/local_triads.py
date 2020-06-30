import networkx as nx
import utils
import time
from math import sqrt


def update_T_in(graph, node, community, T_in_prev):
	"""updates the T_in measure if $node is joining the community, given the previous T_in value.

	Args:
		graph ([nx.Graph]): [the given network]
		node ([int]): [the candidate node that is going to join the community in the next step]
		community ([list]): [list of discovered node in with the community index]
		T_in_prev ([int]): [the previous T_in measure of the community without $node]

	Returns:
		[int]: [the updated T_in measure]
	"""
	increment = 0.0
	neighbors = list(graph.neighbors(node))

	for i in range(0, len(neighbors)):
		for j in range(i+1, len(neighbors)):
			if graph.has_edge(neighbors[i], neighbors[j]):
				if neighbors[i] in community and neighbors[j] in community:
					increment += 1

	return T_in_prev + increment


def update_T_ex(graph, node, community, T_ex_prev):
	"""updates the T_ex measure if $node is joining the community, given the previous T_ex value.

	Args:
		graph ([nx.Graph]): [the given network]
		node ([int]): [the candidate node that is going to join the community in the next step]
		community ([list]): [list of discovered node in with the community index]
		T_in_prev ([int]): [the previous T_ex measure of the community without $node]

	Returns:
		[int]: [the updated T_ex measure]
	"""
	increment = 0.0
	decrement = 0.0
	neighbors = list(graph.neighbors(node))

	for i in range(0, len(neighbors)):
		for j in range(i+1, len(neighbors)):
			if graph.has_edge(neighbors[i], neighbors[j]):
				if not neighbors[i] in community and not neighbors[j] in community:
					increment += 1
				elif (neighbors[i] in community) != (neighbors[j] in community):
					decrement += 1

	return T_ex_prev + increment - decrement


def calc_T(T_in, T_ex):
	"""calculates T score based on T_in and T_ex measures.

	Args:
		T_in ([int]): [the given T_in measure]
		T_ex ([int]): [the given T_ex measure]

	Returns:
		[int]: [T score]
	"""
	if T_in > T_ex:
		return T_in * (T_in - T_ex)
	return 0.0


def find_the_best_neighbor(graph, community, shell_set, T_in_prev, T_ex_prev):
	"""computes T score for each node in the shell set, and returns the one with the highest T score.

	Args:
		graph ([nx.Graph]): [the given network]
		community ([list]): [the list of discovered nodes with the same community index]
		shell_set ([list]): [list of direct neighbors of the discovered community]
		T_in_prev ([int]): [the previous T_in measure of the community without any new best candidate joining]
		T_ex_prev ([int]): [the previous T_ex measure of the community without any new best candidate joining]

	Returns:
		[tuple]: [the best node among candidates, its T_in score, and T_ex score]
	"""
	dict_of_T_in = dict()
	dict_of_T_ex = dict()
	dict_of_T = dict()

	for node in shell_set:
		dict_of_T_in[node] = update_T_in(graph, node, community, T_in_prev)
		dict_of_T_ex[node] = update_T_ex(graph, node, community, T_ex_prev)
		dict_of_T[node] = calc_T(dict_of_T_in[node], dict_of_T_ex[node])

	best_node = shell_set[0]
	for node in shell_set[1:]:
		if dict_of_T[node] > dict_of_T[best_node]:
			best_node = node
		elif dict_of_T[node] == dict_of_T[best_node] and dict_of_T_ex[node] < dict_of_T_ex[best_node]:
			best_node = node

	return best_node, dict_of_T_in[best_node], dict_of_T_ex[best_node]


def find_high_degree_neighbor(graph, node):
	"""finds a neighbor of the given node with the highest degree.

	Args:
		graph ([nx.Graph]): [the given network]
		node ([int]): [the given node to find its best neighbor]

	Returns:
		[int]: [a neighbor of the given node with the highest degree]
	"""
	neighbors = list(graph.neighbors(node)) + [node]
	best_neigh = neighbors[0]
	for i in range(1, len(neighbors)):
		if graph.degree[neighbors[i]] > graph.degree[best_neigh]:
			best_neigh = neighbors[i]

	return best_neigh


def community_search(graph, intended_node, start_with_given_node):
	"""	searches for all nodes in the network that should be in the same community as the $intended_node.

	Args:
		graph ([nx.Graph]): [the given network]
		intended_node ([int]): [a node in the network that starts the community expansion]
		start_with_given_node ([bool]): [if Ture: start expansion with the given node, if False: start with a node of highest degree]

	Returns:
		[list]: [the discovered community]
	"""
	initial_node = intended_node
	if not start_with_given_node:
		initial_node = find_high_degree_neighbor(graph, intended_node)
	
	community = [initial_node]

	shell_set = list(graph.neighbors(initial_node))
	if initial_node in shell_set:
		shell_set.remove(initial_node)

	T_in = update_T_in(graph, initial_node, community=[], T_in_prev=0.0)
	T_ex = update_T_ex(graph, initial_node, community=[], T_ex_prev=0.0)
	curr_T = calc_T(T_in, T_ex)

	while len(community) < graph.number_of_nodes() and shell_set != []:

		# find the neighbor with the highest improvement to the T score
		best_node, new_T_in, new_T_ex = find_the_best_neighbor(graph, community, shell_set, T_in, T_ex)
		
		# calculate the new T score
		new_T = calc_T(new_T_in, new_T_ex)

		if new_T >= curr_T:
			T_in, T_ex, curr_T = new_T_in, new_T_ex, new_T
			community.append(best_node)
	
			shell_set.extend(list(set(graph.neighbors(best_node)) - set(community)))
			shell_set = list(set(shell_set))
			shell_set.remove(best_node)

		else:
			break

	if not intended_node in community:
		community = []

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
		initial_node = utils.find_node_highest_degree(graph, all_communities)
		community = [initial_node]
		nodes_discovered += 1
		
		shell_set = list(graph.neighbors(initial_node))
		if initial_node in shell_set:
			shell_set.remove(initial_node)

		for node in [dummy_node for community in all_communities for dummy_node in community]:
			if node in shell_set:
				shell_set.remove(node)

		T_in = update_T_in(graph, initial_node, community=[], T_in_prev=0.0)
		T_ex = update_T_ex(graph, initial_node, community=[], T_ex_prev=0.0)
		curr_T = calc_T(T_in, T_ex)

		while len(community) < graph.number_of_nodes() and shell_set != []:

			# find the neighbor with the highest improvement to the T score
			best_node, new_T_in, new_T_ex = find_the_best_neighbor(graph, community, shell_set, T_in, T_ex)

			# calculate the new T score
			new_T = calc_T(new_T_in, new_T_ex)

			if new_T >= curr_T:
				T_in, T_ex, curr_T = new_T_in, new_T_ex, new_T
				community.append(best_node)
				nodes_discovered += 1

				new_neighbors = list(set(graph.neighbors(best_node)) - set(community))
				for node in [dummy_node for community in all_communities for dummy_node in community]:
					if node in new_neighbors:
						new_neighbors.remove(node)

				shell_set.extend(new_neighbors)
				shell_set = list(set(shell_set))
				shell_set.remove(best_node)
			
			else:
				break

		all_communities.append(sorted(community))

	# all_communities = utils.amend_partition(graph, all_communities)
	report = utils.report_performance(graph, all_communities, ground_truth_file_address)
	print(report, '\tCommunity Detection Task is done in %.4f seconds.' %(time.time() - start_time))
	return all_communities


def community_search_for_all_nodes(graph, ground_truth_file_address, start_with_given_node=True):
	"""do community search for all nodes separately, calculate accuracy measures then report the Avg. and SD.

	Args:
		graph ([nx.Graph]): [the given network]
		ground_truth_file_address ([str]): [filename of the ground-truth information of communities]
		start_with_given_node (bool, optional): [if Ture: start expansion with the given node, if False: start with a node of highest degree]. Defaults to False.
	"""
	start_time = time.time()
	performance_info = dict()
	ground_truth_com2nodes = utils.read_ground_truth(ground_truth_file_address)

	for e, node in enumerate(graph.nodes()):
		performance_info[node] = {'degree': graph.degree[node], 'precision': 0.0, 'recall': 0.0, 'f1-score': 0.0}
		community = community_search(graph, node, start_with_given_node)
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
