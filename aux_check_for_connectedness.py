import networkx as nx
import utils


def read_communities(filename_address):
	communities = dict()
	with open(filename_address, 'r') as file:
		lines = file.readlines()
		for line in lines:
			line = line.split()
			if not int(line[1]) in communities:
				communities[int(line[1])] = list()
			communities[int(line[1])].append(int(line[0]))

	return communities


def _main():
	args = utils.create_argument_parser()
	graph = utils.load_graph(args.dataset, args.w)
	communities = read_communities(args.output)

	num_of_disconnected_coms = 0
	for com_index, nodes in communities.items():
		small_graph = graph.subgraph(nodes)
		if not nx.is_connected(small_graph):
			num_of_disconnected_coms += 1

	print('\nPercent of disconnected communities = %.3f' %(num_of_disconnected_coms / len(communities)))
	print(num_of_disconnected_coms, '\t', len(communities))


if __name__ == "__main__":
	_main()