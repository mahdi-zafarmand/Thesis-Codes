from networkx.algorithms.community import greedy_modularity_communities
import utils
import time
from measures import modularity
from measures import NMI


def create_partition(communities):
	node2com = dict()
	for com_index, nodes in communities.items():
		for node in nodes:
			node2com[node] = com_index
	return node2com


def main():
	start_time = time.time()

	args = utils.create_argument_parser()
	graph = utils.load_graph(args.dataset, args.w)
	c = greedy_modularity_communities(graph)
	communities = dict()

	for i in range(len(c)):
		communities[i] = list(c[i])

	partition = create_partition(communities)
	utils.print_comm_info_to_display(communities)
	# utils.write_comm_info_to_file(partition)

	print('modularity_value =', modularity(graph, communities))
	print('NMI =', NMI(args.output, partition))

	finish_time = time.time()
	print('\nDone in %.4f seconds.' %(finish_time - start_time))




if __name__ == "__main__":
	main()

