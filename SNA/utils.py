from random import choice


def remove_self_loops(graph):
	for node in graph.nodes():
		if graph.has_edge(node, node):
			graph.remove_edge(node, node)


def find_random_node_in_graph(graph, sets_of_nodes_to_ignore):
	nodes = list(graph.nodes())
	for a_set_of_nodes_to_ignore in sets_of_nodes_to_ignore:
		for node in a_set_of_nodes_to_ignore:
			nodes.remove(node)
	return choice(nodes)


def find_highest_degree_node_in_graph(graph, sets_of_nodes_to_ignore):
	nodes = list(graph.nodes())
	for a_set_of_nodes_to_ignore in sets_of_nodes_to_ignore:
		for node in a_set_of_nodes_to_ignore:
			nodes.remove(node)

	node_with_highest_degree = nodes[0]
	for node_index in range(1, len(nodes)):
		if graph.degree[nodes[node_index]] > graph.degree[node_with_highest_degree]:
			node_with_highest_degree = nodes[node_index]

	return node_with_highest_degree
