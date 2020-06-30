import networkx as nx
import utils
import time
from copy import deepcopy
import local_triads
from math import sqrt


_N_STEP = 3
_K_PARAM = 3


def insise_triangles(graph, community):
	count = 0
	if len(community) < 3:
		return count

	for i in range(0, len(community)):
		for j in range(i+1, len(community)):
			for k in range(j+1, len(community)):
				if graph.has_edge(community[i], community[j]):
					if graph.has_edge(community[i], community[k]):
						if graph.has_edge(community[j], community[k]):
							count += 1
	return count


def outside_triangles(graph, community):
	count = 0
	for node in community:
		neighbors = list(set(graph.neighbors(node)) - set(community))
		if len(neighbors) < 2:
			continue

		for i in range(0, len(neighbors)):
			for j in range(i+1, len(neighbors)):
				if graph.has_edge(neighbors[i], neighbors[j]):
						count += 1
	return count


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


def find_next_neighbor(graph, node, depth=2):
	if _N_STEP == 1:
		return node
	
	if depth < _N_STEP:
		the_neighbors = dict()
		for neigh in graph.neighbors(node):
			the_neighbors[neigh] = find_next_neighbor(graph, neigh, depth + 1)

	else:
		the_neighbors = list(graph.neighbors(node))

	return the_neighbors


def read_trajectory(trajectory):
	type_values = type(list(trajectory.values())[0])
	if type_values == dict:
		output = dict()
		for k,v in trajectory.items():
			output[k] = read_trajectory(trajectory[k])
		return output
	elif type_values == list:
		output = list()
		for k, v in trajectory.items():
			for v_i in v:
				output.append((k, v_i))
		return output
	else:
		return list(trajectory.values())

def clear_invalid_states(forward_states, community):
	i = 0
	while i < len(forward_states):
		if len(set(forward_states[i])) != len(forward_states[i]):
			forward_states.pop(i)
		elif len(set(forward_states[i]).intersection(community)) > 0:
			forward_states.pop(i)
		else:
			i += 1
	return forward_states


def find_forward_states(graph, community):
	trajectory = dict()
	neighbors = list()
	for node in community:
		neighbors.extend(graph.neighbors(node))
	neighbors = list(set(neighbors))

	for neigh in neighbors:
		trajectory[neigh] = find_next_neighbor(graph, neigh)

	states = read_trajectory(trajectory)
	while type(states) != list:
		states_copy = deepcopy(states)
		states = read_trajectory(states_copy)

	forward_states = list()
	if type(states[0]) == int:
		forward_states = [[x] for x in states]
	else:
		temp = list()
		for state in states:
			while type(state[1]) != int:
				temp.append(state[0])
				state = state[1]
			temp.append(state[0])
			temp.append(state[1])
			forward_states.append(temp[:])
			temp.clear()

	return clear_invalid_states(forward_states, community)


def compute_forward_scores(graph, T_in, T_ex, community, forward_states):
	forward_scores = dict()
	internal_scores = dict()
	external_scores = dict()

	for state in forward_states:
		T_in_state = T_in
		T_ex_state = T_ex
		for next_node in state:
			T_in_state = update_T_in(graph, next_node, community, T_in_state)
			T_ex_state = update_T_ex(graph, next_node, community, T_ex_state)
			community.append(next_node)

		for i in range(_N_STEP):
			community.pop()

		internal_scores[tuple(state)] = T_in_state
		external_scores[tuple(state)] = T_ex_state
		forward_scores[tuple(state)] = calc_T(T_in_state, T_ex_state)

	return forward_scores, internal_scores, external_scores


def find_best_state(forward_scores, external_scores):
	states = list(forward_scores.keys())

	best_state = states[0]
	for i in range(1, len(states)):
		if forward_scores[states[i]] > forward_scores[states[i]]:
			best_state = states[i]
		elif forward_scores[states[i]] == forward_scores[states[i]]:
			if external_scores[states[i]] < external_scores[best_state]:
				best_state = states[i]

	return best_state


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
	initial_node = intended_node
	if not start_with_given_node:
		initial_node = find_high_degree_neighbor(graph, intended_node)

	community = list()

	community.append(initial_node)
	T_in = update_T_in(graph, initial_node, [], 0.0)
	T_ex = update_T_ex(graph, initial_node, [], 0.0)
	curr_T = calc_T(T_in, T_ex)
	k_param_counter = 0
	should_stop = False

	while len(community) < graph.number_of_nodes() and k_param_counter < _K_PARAM:
		k_param_counter += 1

		forward_states = find_forward_states(graph, community)
		forward_scores, internal_scores, external_scores = compute_forward_scores(graph, T_in, T_ex, community, forward_states)
		if len(forward_scores) == 0:
			break
		best_state = find_best_state(forward_scores, external_scores)
		best_next_node = best_state[0]

		new_T_in = update_T_in(graph, best_next_node, community, T_in)
		new_T_ex = update_T_ex(graph, best_next_node, community, T_ex)
		new_T = calc_T(new_T_in, new_T_ex)

		if new_T >= curr_T:
			T_in, T_ex, curr_T = new_T_in, new_T_ex, new_T
			community.append(best_next_node)
		else:
			should_stop = True
			break

	shell_set = set()
	for node in community:
		shell_set.update(set(graph.neighbors(node)))
	shell_set = list(shell_set - set(community))

	while len(community) < graph.number_of_nodes() and shell_set != [] and should_stop == False:
		best_node, new_T_in, new_T_ex = local_triads.find_the_best_neighbor(graph, community, shell_set, T_in, T_ex)		
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
	start_time = time.time()
	graph_copy = deepcopy(graph)
	all_communities = list()

	while graph_copy.number_of_nodes() > 0:
		nodes = list(graph_copy.nodes())
		intended_node = nodes[0]
		for i in range(1, len(nodes)):
			if graph_copy.degree[nodes[i]] > graph_copy.degree[intended_node]:
				intended_node = nodes[i]

		community = community_search(graph_copy, intended_node, True)
		all_communities.append(sorted(community))
		graph_copy.remove_nodes_from(community)

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
