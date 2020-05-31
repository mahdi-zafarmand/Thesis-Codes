from networkx.algorithms.community import greedy_modularity_communities
import utils
from copy import deepcopy
import time
from measures import modularity
from measures import NMI


def create_partition(communities):
	node2com = dict()
	for com_index, nodes in communities.items():
		for node in nodes:
			node2com[node] = com_index
	return node2com


def preprocess(graph_base, weight="weight"):
	nodes = list(graph_base.nodes())
	nb_vertices = len(nodes)
	mutuals, max_mutuals, total_mutuals = get_mutuals(graph_base)

	edges_done = list()

	for i in range(nb_vertices):
		neighbors = list(graph_base.neighbors(nodes[i]))
		max_mutual_node = max_mutuals.get(nodes[i])

		for neigh in neighbors:
			if (neigh, nodes[i]) in edges_done:
				continue
			max_mutual_neigh = max_mutuals.get(neigh)
			w = mutuals.get(nodes[i]).get(neigh)
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
			w = (w1 + w2) / 2
			graph_base.add_edge(nodes[i], neigh, weight= w)
			edges_done.append((nodes[i], neigh))


def get_mutuals(graph_base):
	mutuals = {}		# mutual[n1][n2]	= number of mutual neighbors between nodes n1 and n2
	max_mutuals = {}	# max_mutuals[n]	= maximum number of mutuals between node n and any other node = max(mutuals[n][k])
	total_mutuals = {}	# total_mutuals[n]	= number of total mutual of node n with any other node = signma(mutuals[n][k])
	nodes = list(graph_base.nodes())
	nb_vertices = len(nodes)

	# initializing the return values
	for i in range(nb_vertices):
		mutuals[nodes[i]] = {}
		max_mutuals[nodes[i]] = -1
		total_mutuals[nodes[i]] = 0.0

	for i in range(nb_vertices):
		neighbors = list(graph_base.neighbors(nodes[i]))
		for neigh in neighbors:
			if neigh in mutuals[nodes[i]]:
				continue
			cur_mutual = shared_neighbors_cnt(graph_base, nodes[i], neigh)
			mutuals[nodes[i]][neigh] = cur_mutual
			mutuals[neigh][nodes[i]] = cur_mutual

			total_mutuals[nodes[i]] = total_mutuals[nodes[i]] + cur_mutual
			total_mutuals[neigh] = total_mutuals[neigh] + cur_mutual

			if cur_mutual > max_mutuals[nodes[i]]:
				max_mutuals[nodes[i]] = cur_mutual
			if cur_mutual > max_mutuals[neigh]:
				max_mutuals[neigh] = cur_mutual

		# to avoid the double-counting
		total_mutuals[nodes[i]] = total_mutuals[nodes[i]]/2.0

	return mutuals, max_mutuals, total_mutuals


def shared_neighbors_cnt(graph_base, u, v):
	shared = 0
	if(graph_base.degree(u) > graph_base.degree(v)):
		tmp = u
		u = v
		v = tmp
	neighbors_u = graph_base[u]
	neighbors_v = graph_base[v]
	for n1 in neighbors_u:
		if n1 in neighbors_v:
			shared = shared+1
	return shared


def main():
	start_time = time.time()

	args = utils.create_argument_parser()
	graph = utils.load_graph(args.dataset, args.w)
	graph_copy = deepcopy(graph)

	preprocess(graph)

	c = greedy_modularity_communities(graph)
	communities = dict()

	for i in range(len(c)):
		communities[i] = list(c[i])

	partition = create_partition(communities)
	utils.print_comm_info_to_display(communities)
	# utils.write_comm_info_to_file(partition)

	print('modularity_value =', modularity(graph_copy, communities))
	print('NMI =', NMI(args.output, partition))

	finish_time = time.time()
	print('\nDone in %.4f seconds.' %(finish_time - start_time))


if __name__ == "__main__":
	main()

