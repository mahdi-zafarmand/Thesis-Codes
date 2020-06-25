import networkx as nx
import argparse
import utils
import time


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
		w = w1 + w2 - 1.0
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

	return quality / 2.0


def community_search(graph, intended_node):
	mutuals = dict()
	max_mutuals = dict()

	community = list()
	cur_quality = 0.0

	community.append(intended_node)
	assign_local_strength(graph, intended_node, mutuals, max_mutuals)

	while len(community) < graph.number_of_nodes():
		new_candidates = find_new_candidates(graph, community)

		if new_candidates == []:
			break

		for candidate in new_candidates:
			assign_local_strength(graph, candidate, mutuals, max_mutuals)

		new_node = expand_community(graph, community, new_candidates)
		new_quality = compute_total_strength(graph, community, new_node)

		if new_quality - cur_quality < __MIN:
			break

		community.append(new_node)
		cur_quality = new_quality

	return community


def amend_by_dangles(graph, community):
	dangles = list()
	for node in community:
		for neigh in graph.neighbors(node):
			if graph.degree[neigh] == 1:
				dangles.append(neigh)

	return sorted(community + dangles)


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

	for e, node in enumerate(graph.nodes()):
		community = community_search(graph, node)
		community = amend_by_dangles(graph, community)
		communities[node] = community

	ground_truth = read_ground_truth(args.output)
	precision, recall, f1_score = calc_accuracy(communities, ground_truth)

	print('precision =', precision, '\trecall =', recall, '\tf1-score =', f1_score)

	finish_time = time.time()
	print('\nDone in %.4f seconds.' %(finish_time - start_time))


if __name__ == "__main__":
	_main()
