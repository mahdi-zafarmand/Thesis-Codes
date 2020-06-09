import networkx as nx
from measures import modularity
import utils


def _main():
	args = utils.create_argument_parser()
	graph = utils.load_graph(args.dataset, args.w, self_loop=True)

	nodes = list(graph.nodes()) 
	print('num of nodes =', len(nodes))
	print('num of edges =', graph.number_of_edges())

	communities = dict()
	communities_file = args.output
	with open(communities_file, 'r') as file:
		lines = file.readlines()
		for line in lines:
			line = line.split()
			if not int(line[0]) in nodes:
				continue
			if not int(line[1]) in communities:
				communities[int(line[1])] = list()
			communities[int(line[1])].append(int(line[0]))

	print('modularity_value =', modularity(graph, communities))


if __name__ == "__main__":
	_main()