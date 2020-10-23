from copy import deepcopy
import utils


class CommunitySearcher:
	# the high-level class to search for the community of a given node in a social network.
	minimum_improvement = 0.000001     # if an improvement is less than minimum, the process stops considering stability

	def __init__(self, name, graph):
		# initializes the object, not all attributes are useful for every algorithm
		self.name = name
		self.graph = deepcopy(graph)
		self.nodes_to_be_ignored = set()
		self.starting_node = None
		self.community = []
		self.boundary = set()
		self.shell = set()
		self.remove_self_loops()

	def reset(self):
		# clears auxiliary sets in order to use the same CommunitySearcher object to find next communities.
		self.community.clear()
		self.boundary.clear()
		self.shell.clear()

	def fill_ignored_nodes(self, nodes):
		# to populate the nodes that should be ignored when candidates are being picked.
		self.nodes_to_be_ignored.update(set(nodes))

	def empty_ignored_nodes(self):
		# should be called when we need community detection with overlap.
		self.nodes_to_be_ignored.clear()

	def remove_self_loops(self):
		# algorithms tend to work better if there is no self-loop in the given graph, so we call this method at first.
		utils.remove_self_loops(self.graph)

	def set_start_node(self, start_node):
		# check the validity of the given start_node, then puts it in the community and initialize the shell set.
		if start_node in self.graph.nodes():
			self.starting_node = start_node
			self.community.append(start_node)
			self.find_initial_shell_set()
		else:
			print('Invalid starting node! Try with another one.')
			exit(-1)

	def find_initial_shell_set(self):
		# constructs the initial shell set, which contains the candidates for the next node to join the community.
		self.shell = set(self.graph.neighbors(self.starting_node))
		for node in self.nodes_to_be_ignored:
			if node in self.shell:
				self.shell.remove(node)

	def community_search(self, start_node, with_amend=False):
		# has different implementation for different algorithms, using "polymorphism" and "method overriding".
		pass

	def compute_modularity(self, neighbor_node):
		# has different implementation for different algorithms, using "polymorphism" and "method overriding".
		pass

	def update_boundary(self, new_node):
		# after a new_node expands the community, boundary set should be updated by adding and removing some nodes.
		neighbors_of_new_node = self.graph.neighbors(new_node)
		should_be_boundary = False
		for neighbor in neighbors_of_new_node:
			if (neighbor in self.community) is False:
				should_be_boundary = True
				break
		if should_be_boundary:
			self.boundary.add(new_node)

		# if the only neighbor of a node that was out of community was new_node, then it should leave the boundary.
		possibles_leaving_nodes = [node for node in neighbors_of_new_node if node in self.boundary]
		for node in possibles_leaving_nodes:
			should_leave_boundary = True
			for neighbor in self.graph.neighbors(node):
				if (neighbor in self.community) is False:
					should_leave_boundary = False
					break
			if should_leave_boundary:
				self.boundary.remove(node)

	def update_shell(self, new_node):
		# after a new_node expands the community, the shell set should be updated.
		self.shell.update(self.graph.neighbors(new_node))
		for node in self.community:
			if node in self.shell:
				self.shell.remove(node)

		for node in self.nodes_to_be_ignored:
			if node in self.shell:
				self.shell.remove(node)

	def show_community(self):
		print("There are " + str(len(self.community)) + " nodes in the community.")
		print("C = {" + str(self.community) + "}")
