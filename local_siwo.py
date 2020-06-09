import networkx as nx
import argparse
import utils
import time


__MIN = 0.0000001


def create_argument_parser():
	parser = argparse.ArgumentParser(description='The SIWO algorithm', prog="SIWO")
	parser.add_argument('dataset', help='Path to the network files')
	parser.add_argument('output', help='output name')
	parser.add_argument('support', help='support value')
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
		w = (w1 + w2)/2.0
		graph.add_edge(node, neigh, strength=w)


def find_new_candidates(graph, community):
	candidates = set()
	for node in community:
		for neigh in graph.neighbors(node):
			candidates.add(neigh)

	return list(candidates - set(community))


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
	all_node = community + [new_node]
	quality = 0.0
	for node in all_node:
		for neigh in graph.neighbors(node):
			if neigh in all_node:
				quality += graph[node][neigh].get('strength', 0.0)

	return quality / 2


def community_search_one_round(graph, initial_node):
	cur_quality = 0.0	
	mutuals = {}
	max_mutuals = {}

	community = [initial_node]
	assign_local_strength(graph, initial_node, mutuals, max_mutuals)

	while True:
		new_candidates = find_new_candidates(graph, community)
		new_node = expand_community(graph, community, new_candidates)
		new_quality = compute_total_strength(graph, community, new_node)
		
		if new_quality - cur_quality < __MIN:
			break

		community.append(new_node)
		cur_quality = new_quality

	return sorted(community), cur_quality


# def community_search_multiple(graph, intended_node, min_support):
# 	communities = dict()
# 	community, quality_main = community_search_one_round(graph, intended_node)
# 	communities[intended_node] = (community, quality_main)

# 	for i in range(1, len(community)):
# 		comm, quality_other = community_search_one_round(graph, community[i])
# 		communities[community[i]] = (comm, quality_other)

# 	all_candidates = dict()
# 	for com, (nodes, q) in communities.items():
# 		for node in nodes:
# 			all_candidates[node] = all_candidates.get(node, 0) + 1

# 	final_answer = list()
# 	for candidate, num_repeat in all_candidates.items():
# 		if num_repeat > min_support * len(communities):
# 			final_answer.append(candidate)

# 	return final_answer


def community_search_multiple(graph, intended_node, min_support):
	communities = dict()
	community, quality_main = community_search_one_round(graph, intended_node)
	communities[intended_node] = (community, quality_main)

	neighbors = graph.neighbors(intended_node)

	for neigh in neighbors:
		comm, quality_other = community_search_one_round(graph, neigh)
		communities[neigh] = (comm, quality_other)

	# for k, v in communities.items():
	# 	print(k, sorted(v[0]), v[1], len(v[0]), v[1] / len(v[0]))

	all_candidates = dict()
	for com, (nodes, q) in communities.items():
		for node in nodes:
			all_candidates[node] = all_candidates.get(node, 0) + 1

	final_answer = list()
	for candidate, num_repeat in all_candidates.items():
		if num_repeat > min_support * len(communities):
			final_answer.append(candidate)

	return final_answer



def amend_by_dangles(graph, community):
	dangles = list()
	for node in community:
		for neigh in graph.neighbors(node):
			if graph.degree[neigh] == 1:
				dangles.append(neigh)

	return sorted(community + dangles)


def _main():
	start_time = time.time()

	args = create_argument_parser()
	graph = utils.load_graph(args.dataset, False)

	intended_node = int(args.output)
	min_support = float(args.support)

	community = community_search_one_round(graph, intended_node)[0]
	# community = community_search_multiple(graph, intended_node, min_support)
	# community = amend_by_dangles(graph, community)
	print('community =', community, len(community))

	finish_time = time.time()
	print('\nDone in %.4f seconds.' %(finish_time - start_time))


if __name__ == "__main__":
	_main()