import networkx as nx
import argparse
import utils
import time
from measures import modularity
from measures import NMI


__MIN = 0.0000001


def create_argument_parser():
	parser = argparse.ArgumentParser(description='The SIWO algorithm', prog="SIWO")
	parser.add_argument('dataset', help='Path to the network files')
	parser.add_argument('output', help='output name')
	return parser.parse_args()


def shared_neighbors_cnt(graph_base, u, v):
	shared = 0
	if(graph_base.degree(u) > graph_base.degree(v)):
		u, v = v, u
	neighbors_u = graph_base[u]
	neighbors_v = graph_base[v]
	for n1 in neighbors_u:
		if n1 in neighbors_v:
			shared = shared+1
	return shared


def update_mutual_dicts(graph, node, mutuals, max_mutuals):
	neighbors = list(graph.neighbors(node))
	if not node in mutuals:
		mutuals[node] = {}
		max_mutuals[node] = -1

	for neigh in neighbors:
		if neigh in mutuals[node]:
			continue
		if not neigh in mutuals:
			mutuals[neigh] = {}
			max_mutuals[neigh] = -1

		cur_mutual = shared_neighbors_cnt(graph, node, neigh)
		mutuals[node][neigh] = cur_mutual
		mutuals[neigh][node] = cur_mutual

		if cur_mutual > max_mutuals[node]:
			max_mutuals[node] = cur_mutual
		if cur_mutual > max_mutuals[neigh]:
			max_mutuals[neigh] = cur_mutual


def assign_local_strength(graph, node, mutuals, max_mutuals):
	update_mutual_dicts(graph, node, mutuals, max_mutuals)
	max_mutual_node = max_mutuals.get(node)

	for neigh in graph.neighbors(node):
		max_mutual_neigh = max_mutuals.get(neigh)
		w = mutuals.get(node).get(neigh)
		w1 = 0.0
		w2 = 0.0
		try:
			w1 = w / max_mutual_node
		except:
			pass
		try:
			w2 = w / max_mutual_neigh
		except:
			pass
		w = (w1 + w2) - 1.0
		graph.add_edge(node, neigh, strength=w)


def find_new_candidates(graph, community):
	candidates = set()
	for node in community:
		for neigh in graph.neighbors(node):
			candidates.add(neigh)

	return list(candidates - set(community))


def clear_former_nodes(candidates, communities):
	if len(communities) > 0:
		for comm in communities:
			for node in comm:
				if node in candidates:
					candidates.remove(node)
	return candidates


def expand_community(graph, community, candidates):
	sum_strength = dict()
	for candidate in candidates:
		sum_strength[candidate] = 0.0
		for neigh in graph.neighbors(candidate):
			if not neigh in community:
				continue
			sum_strength[candidate] += graph[candidate][neigh].get('strength', 0.0)

	best_candidate = candidates[0]
	for candidate in candidates:
		if sum_strength[candidate] > sum_strength[best_candidate]:
			best_candidate = candidate

	return best_candidate

def compute_total_strength(graph, community, new_node):
	new_community = community + [new_node]
	quality = 0.0
	for node in new_community:
		for neigh in graph.neighbors(node):
			if neigh in new_community:
				quality += graph[node][neigh].get('strength', 0.0)

	return quality / 2.0


def find_node_highest_degree(graph, communities):
	nodes = list(graph.nodes())

	for comm in communities:
		for node in comm:
			if node in nodes:
				nodes.remove(node)

	highest_deg_node = nodes[0]
	for i in range(1, len(nodes)):
		if graph.degree[nodes[i]] > graph.degree[highest_deg_node]:
			highest_deg_node = nodes[i]

	return highest_deg_node



def community_search(graph):
	communities = list()
	mutuals = dict()
	max_mutuals = dict()

	nodes_discovered = 0

	while nodes_discovered < graph.number_of_nodes():
		cur_quality = 0.0	
		community = list()
		initial_node = find_node_highest_degree(graph, communities)			# among undiscovered nodes
		community.append(initial_node)
		nodes_discovered += 1
		assign_local_strength(graph, initial_node, mutuals, max_mutuals)	# of the edges that are connected to the initial node

		while True:
			new_candidates = find_new_candidates(graph, community)			# the direct neighbors of the community
			new_candidates = clear_former_nodes(new_candidates, communities)# remove all nodes that have been discovered
			
			if new_candidates == []:
				break

			for candidate in new_candidates:								# to the edges that connects the community to its neighborhood
				assign_local_strength(graph, candidate, mutuals, max_mutuals)

			new_node = expand_community(graph, community, new_candidates)	# find the best candidate to merge into the community
			new_quality = compute_total_strength(graph, community, new_node)# compute new quality if the best candidate joins

			if new_quality - cur_quality < __MIN:							# if improvement is less than threshold (zero), stop!
				break				

			community.append(new_node)										# adds the best candidates among all neighbors to the community
			nodes_discovered += 1
			cur_quality = new_quality

		communities.append(sorted(community))								# sort is for better representation only!
	
	return communities


def amend_singleton_communities(graph, communities):
	singleton_coms_indices = [i for i in range(len(communities)) if len(communities[i]) == 1]
	best_neigh_coms_indices = dict()

	for i in singleton_coms_indices:
		neighbors = list(graph.neighbors(communities[i][0]))
		neigh_comms = dict()

		for j in range(len(communities)):
			neigh_comms[j] = len([x for x in neighbors if x in communities[j]])

		best_com = 0
		for j in range(1, len(communities)):
			if neigh_comms[j] > neigh_comms[best_com]:
				best_com = j

		best_neigh_coms_indices[i] = best_com

	for old_com_index, new_com_index in best_neigh_coms_indices.items():
		communities[new_com_index].append(communities[old_com_index][0])

	best_neigh_coms_indices = list(best_neigh_coms_indices)
	best_neigh_coms_indices = sorted(best_neigh_coms_indices, reverse=True)
	for i in range(len(best_neigh_coms_indices)):
		communities.pop(best_neigh_coms_indices[i])

	return communities


def _main():
	start_time = time.time()

	args = create_argument_parser()
	graph = utils.load_graph(args.dataset, False)

	communities = community_search(graph)
	communities = amend_singleton_communities(graph, communities)
	com_dict = {}

	for i in range(len(communities)):
		com_dict[i] = communities[i]

	utils.print_comm_info_to_display(com_dict)
	print('modularity_value =', modularity(graph, com_dict))

	com_dict2 = {}
	for k, v in com_dict.items():
		for node in v:
			com_dict2[node] = k
			
	print('NMI =', NMI(args.output, com_dict2))

	finish_time = time.time()
	print('\nDone in %.4f seconds.' %(finish_time - start_time))


if __name__ == "__main__":
	_main()