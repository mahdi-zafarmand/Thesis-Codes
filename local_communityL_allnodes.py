import networkx as nx
import utils
import time


def compute_L(graph, community, boundary, new_node):
	new_boundary = boundary + [new_node]
	new_community = community + [new_node]
	L_in = 0.0
	L_ex = 0.0

	for node in new_community:
		for neigh in graph.neighbors(node):
			if neigh in new_community:
				L_in += 1
	L_in /= len(new_community)

	for node in new_boundary:
		for neigh in graph.neighbors(node):
			if not neigh in new_boundary and not neigh in community:
				L_ex += 1
	L_ex /= len(new_boundary)

	return L_in / L_ex	


def comupte_improvements_for_neighborhood(graph, community, boundary, neighbors, L_measure):
	delta_L = dict()
	for neigh in neighbors:
		delta_L[neigh] = compute_L(graph, community, boundary, neigh) - L_measure

	return delta_L


def find_the_best_neighbor(delta_L):
	neighbors = list(delta_L.keys())
	best_neigh = neighbors[0]
	for i in range(1, len(neighbors)):
		if delta_L[neighbors[i]] > delta_L[best_neigh]:
			best_neigh = neighbors[i]

	return best_neigh


def update_neighbors(graph, community, the_neighbors, best_neigh):
	neighbors_of_best_neigh = list(set(graph.neighbors(best_neigh)) - set(community))
	if neighbors_of_best_neigh != []:
		the_neighbors.extend(neighbors_of_best_neigh)
		the_neighbors = list(set(the_neighbors))
	the_neighbors.remove(best_neigh)
	return the_neighbors


def update_boundaries(graph, community, boundary, best_neigh):
	# add best_neigh to the boundary if necessary 
	# (if the best_neigh is only connected to the core, it won't be added to boundaries)
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
	community = list()
	boundary = list()
	neighbors = list()
	L_measure = 0.0

	community.append(intended_node)
	boundary.append(intended_node)
	neighbors = list(graph.neighbors(intended_node))

	while len(community) < graph.number_of_nodes():

		# compute the improvement caused by any node in the neighborhood
		delta_L = comupte_improvements_for_neighborhood(graph, community, boundary, neighbors, L_measure)

		# find the neighbor with the highest improvement to the modularity R
		best_neigh = find_the_best_neighbor(delta_L)

		# continue expanding the community only if modularity R increases
		if delta_L[best_neigh] < 0.0:
			break

		# update modularity R and add the best neighbor to the community
		L_measure += delta_L[best_neigh]
		community.append(best_neigh)

		# update the neighborhood list
		neighbors = update_neighbors(graph, community, neighbors, best_neigh)

		# update the boundary list
		boundary = update_boundaries(graph, community, boundary, best_neigh)

	return sorted(community)



def read_ground_truth(file_address):
	ground_truth = dict()
	with open(file_address, 'r') as file:
		lines = file.readlines()
		for line in lines:
			line = line.split()
			ground_truth[int(line[0])] = int(line[1])

	communities = dict()
	for node1 in ground_truth:
		new_set = set()
		for node2 in ground_truth:
			if ground_truth[node1] == ground_truth[node2]:
				new_set.add(node2)
		communities[node1] = list(new_set)

	return communities


def calc_accuracy(communities, ground_truth):
	precision = dict()
	recall = dict()
	f1_score = dict()

	for key in communities.keys() & ground_truth.keys():
		precision[key] = len(set(communities[key]) & set(ground_truth[key]))

	for key in precision:
		try:
			recall[key] = precision[key] / len(ground_truth[key])
		except:
			recall[key] = 0.0
		try:
			precision[key] = precision[key] / len(communities[key])
		except:
			precision[key] = 0.0
		try:
			f1_score[key] = 2 * (recall[key] * precision[key]) / (recall[key] + precision[key])
		except:
			f1_score[key] = 0.0


	avg_precision = 0.0
	avg_recall = 0.0
	avg_f1score = 0.0

	for key in precision:
		avg_precision += precision[key]
		avg_recall += recall[key]
		avg_f1score += f1_score[key]

	return avg_precision / len(communities) , avg_recall / len(communities), avg_f1score / len(communities)


def _main():
	start_time = time.time()

	args = utils.create_argument_parser()
	graph = utils.load_graph(args.dataset, args.w)

	communities = dict()

	print('\n\n')
	for e, node in enumerate(sorted(graph.nodes())):
		community = community_search(graph, node)
		communities[node] = community
		# print('  Node =', node, ' : degree =', graph.degree[node], '->\t', community, '\t', len(community))

	ground_truth = read_ground_truth(args.output)
	precision, recall, f1_score = calc_accuracy(communities, ground_truth)

	print('precision =', precision, '\trecall =', recall, '\tf1-score =', f1_score)

	finish_time = time.time()
	print('\nDone in %.4f seconds.' %(finish_time - start_time))


if __name__ == "__main__":
	_main()