import utils
import time
import networkx as nx
import random
import os.path
from operator import eq


def calc_T_in(graph, node, community, Tin_prev):
	increment = 0.0
	neighbors = list(graph.neighbors(node))

	for neigh1 in neighbors:
		for neigh2 in neighbors:
			if neigh1 != neigh2 and graph.has_edge(neigh1, neigh2):
				if neigh1 in community and neigh2 in community:
					increment += 1

	return Tin_prev + 0.5 * increment


def calc_T_out(graph, node, community, Tout_prev):
	increment = 0.0
	decrement = 0.0
	neighbors = list(graph.neighbors(node))

	for neigh1 in neighbors:
		for neigh2 in neighbors:
			if neigh1 != neigh2 and graph.has_edge(neigh1, neigh2):
				if not neigh1 in community and not neigh2 in community:
					increment += 1
				elif neigh1 in community and not neigh2 in community:
					decrement += 1

	return Tout_prev + 0.5 * increment - decrement


def calc_T(Tin, Tout):
	if Tin > Tout:
		return Tin * (Tin - Tout)
	return 0.0


def find_best_next_node(graph, community, shell_set, Tin_prev, Tout_prev):
	dict_of_Tin = dict()
	dict_of_Tex = dict()
	dict_of_T = dict()

	for node in shell_set:
		T_in_node = calc_T_in(graph, node, community, Tin_prev)
		T_out_node = calc_T_out(graph, node, community, Tout_prev)
		dict_of_Tin[node] = T_in_node
		dict_of_Tex[node] = T_out_node
		dict_of_T[node] = calc_T(T_in_node, T_out_node)

	best_node = shell_set[0]
	for i in range(1, len(shell_set)):
		if dict_of_T[shell_set[i]] > dict_of_T[best_node]:
			best_node = shell_set[i]
		elif dict_of_T[shell_set[i]] == dict_of_T[best_node]:
			if dict_of_Tex[shell_set[i]] < dict_of_Tex[best_node]:
				best_node = shell_set[i]

	return best_node, dict_of_Tin[best_node], dict_of_Tex[best_node]


def find_best_neighbor(graph, node):
	neighbors = list(graph.neighbors(node)) + [node]
	best_neigh = neighbors[0]

	for i in range(1, len(neighbors)):
		if graph.degree[neighbors[i]] > graph.degree[best_neigh]:
			best_neigh = neighbors[i]

	return best_neigh


def community_search(graph, given_node):
	initial_node = find_best_neighbor(graph, given_node)
	community = [initial_node]
	shell_set = list(graph.neighbors(initial_node))
	Tin = calc_T_in(graph, initial_node, [], 0.0)
	Tout = calc_T_out(graph, initial_node, [], 0.0)

	while len(community) < graph.number_of_nodes():
		best_node, new_Tin, new_Tout = find_best_next_node(graph, community, shell_set, Tin, Tout)
		T = calc_T(new_Tin, new_Tout)
		if T >= calc_T(Tin, Tout):
			Tin, Tout = new_Tin, new_Tout
			new_neighbors = list(set(graph.neighbors(best_node)) - set(community))    
			community.append(best_node)
			shell_set.extend(new_neighbors)
			shell_set = list(set(shell_set))
			shell_set.remove(best_node)

			if not shell_set:
				break
		
		else:
			break

	if not given_node in community:
		community = []
		
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


def main():
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
	main()