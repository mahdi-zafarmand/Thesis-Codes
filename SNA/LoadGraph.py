import networkx as nx
import os.path


def load_graph(path, weighted=False, delimiter='\t', self_loop=False):
	graph = nx.Graph()
	if not os.path.isfile(path):
		print("Error: file " + path + " not found!")
		exit(-1)

	with open(path) as file:
		for line in file.readlines():
			w = 1.0
			line = line.split(delimiter)
			v1 = int(line[0])
			v2 = int(line[1])
			graph.add_node(v1)
			graph.add_node(v2)

			if weighted:
				w = float(line[2])

			if (self_loop and v1 == v2) or (v1 != v2):
				graph.add_edge(v1, v2, weight=w)

	return graph


def load_graph_uncertain(path, delimiter='\t', self_loop=False):
	graph = nx.Graph()
	if not os.path.isfile(path):
		print("Error: file " + path + " not found!")
		exit(-1)

	with open(path) as file:
		for line in file.readlines():
			line = line.split(delimiter)
			v1 = int(line[0])
			v2 = int(line[1])
			graph.add_node(v1)
			graph.add_node(v2)

			w = float(line[2])
			p = float(line[3])

			if (self_loop and v1 == v2) or (v1 != v2):
				graph.add_edge(v1, v2, weight=w, prob=p)

	return graph
