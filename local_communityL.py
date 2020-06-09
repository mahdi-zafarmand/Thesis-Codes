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


def _main():
	start_time = time.time()

	args = utils.create_argument_parser()
	graph = utils.load_graph(args.dataset, False)

	intended_node = int(args.output)

	community = community_search(graph, intended_node)
	print('community =', community, len(community))

	finish_time = time.time()
	print('\nDone in %.4f seconds.' %(finish_time - start_time))


if __name__ == "__main__":
	_main()