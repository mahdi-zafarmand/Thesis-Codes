from CommunitySearch import CommunitySearcher
from CommunityDetection import CommunityDetector
import networkx as nx


class LocalSiwoCommunityDiscovery(CommunitySearcher):
	# the class to search for the community of a given node in a social network using Local SIWO algorithm.

	def __init__(self, graph):
		# initializes the object
		super(LocalSiwoCommunityDiscovery, self).__init__('Local SIWO', graph)
		self.dict_common_neighbors = {}
		self.max_common_neighbors = {}
		self.strength_assigned_nodes = set()  # assigning strength values to the edges that are connected to these nodes.

	def reset(self):
		# resets the object to prepare it for another use
		super(LocalSiwoCommunityDiscovery, self).reset()

	def update_dicts_of_common_neighbors_info(self, node):
		# gathers info about the number/max number of common neighbors of the given node and its neighbors
		# (in another word): counts the number of triangles a node participates in and max number of triangles that
		# this node simultaneously participates in with a neighbor node.

		# initializing for a node that has not been visited yet.
		if (node in self.dict_common_neighbors) is False:
			self.dict_common_neighbors[node] = {}
			self.max_common_neighbors[node] = -1

		for neighbor in self.graph.neighbors(node):
			if neighbor in self.dict_common_neighbors[node]:
				continue
			if (neighbor in self.dict_common_neighbors) is False:
				self.dict_common_neighbors[neighbor] = {}
				self.max_common_neighbors[neighbor] = -1

			number_common_neighbors = sum(1 for _ in nx.common_neighbors(self.graph, node, neighbor))
			self.dict_common_neighbors[node][neighbor] = number_common_neighbors
			self.dict_common_neighbors[neighbor][node] = number_common_neighbors

			if number_common_neighbors > self.max_common_neighbors[node]:
				self.max_common_neighbors[node] = number_common_neighbors
			if number_common_neighbors > self.max_common_neighbors[neighbor]:
				self.max_common_neighbors[neighbor] = number_common_neighbors

	def assign_local_strength(self, node):
		# assigns strength to the edges that are connected to the given node, if the node has not been visited before.
		if node in self.strength_assigned_nodes:
			return

		self.update_dicts_of_common_neighbors_info(node)
		max_mutual_node = self.max_common_neighbors.get(node)

		for neighbor in self.graph.neighbors(node):
			max_mutual_neighbor = self.max_common_neighbors.get(neighbor)
			strength = self.dict_common_neighbors.get(node).get(neighbor)
			try:
				s1 = strength / max_mutual_node
			except ZeroDivisionError:
				s1 = 0.0
			try:
				s2 = strength / max_mutual_neighbor
			except ZeroDivisionError:
				s2 = 0.0

			strength = s1 + s2 - 1.0
			self.graph.add_edge(node, neighbor, strength=strength)
		self.strength_assigned_nodes.add(node)

	def find_best_next_node(self, improvements):
		# updates the improvement that can be achieved by merging a node from shell to community, then returns the max.
		new_node = self.community[-1]   # we only update improvements that are affected by adding the last node
		for node in self.shell:
			if (node in improvements) is False:
				improvements[node] = self.graph[node][new_node].get('strength', 0.0)
			elif self.graph.has_edge(node, new_node):
				improvements[node] += self.graph[node][new_node].get('strength', 0.0)
		if new_node in improvements:
			del improvements[new_node]

		best_candidate = None
		best_improvement = -float('inf')
		for candidate in self.shell:
			if improvements[candidate] > best_improvement:
				best_candidate = candidate
				best_improvement = improvements[candidate]

		return best_candidate, best_improvement

	def compute_modularity(self, neighbor_node):
		# not needed.
		pass

	def amend_for_community_search(self):
		# adds any dangling node in the neighborhood of the discovered community.
		neighborhood = set()
		for node in self.community:
			for neighbor in self.graph.neighbors(node):
				neighborhood.add(neighbor)

		dangling_neighbors = [node for node in neighborhood if self.graph.degree[node] == 1]
		self.community = list(set(self.community + dangling_neighbors))

	def community_search(self, start_node, with_amend=False):
		# THE MAIN FUNCTION OF THE CLASS, finds all other nodes that belong to the same community as the start_node does.
		self.set_start_node(start_node)
		self.assign_local_strength(self.starting_node)

		improvements = {}  # key: candidate nodes from the shell set, value: total improved strength after a node joins.
		while len(self.community) < self.graph.number_of_nodes() and self.shell != []:
			for node in self.shell:
				self.assign_local_strength(node)

			new_node, improvement = self.find_best_next_node(improvements)
			if improvement < CommunitySearcher.minimum_improvement:
				break

			self.community.append(new_node)
			self.update_shell(new_node)
			# this algorithm does not use boundary set, so no need to update that.

		if with_amend:
			self.amend_for_community_search()

		return sorted(self.community)   # sort is only for a better representation, can be ignored to boost performance.


class LocalSiwoCommunityDetection(CommunityDetector):
	# the class to detect all communities of a social network by applying Local SIWO algorithm over and over.

	def __init__(self, graph):
		# initialize the object
		super().__init__('Local SIWO', graph)
		self.local_searcher = LocalSiwoCommunityDiscovery(graph)
		self.dict_common_neighbors = {}
		self.max_common_neighbors = {}
		self.strength_assigned_nodes = set()

	def community_detection(self, node_selection_mode, with_amend=False, overlap_enabled=False):
		# THE MAIN FUNCTION OF THE CLASS, finds all communities of the graph.

		while self.number_discovered_nodes < self.graph.number_of_nodes():
			self.select_starting_node(node_selection_mode)
			community = self.local_searcher.community_search(start_node=self.starting_nodes[-1], with_amend=with_amend)
			self.partition.append(community)
			self.local_searcher.reset()
			self.update_number_discovered_nodes()

			if overlap_enabled:
				self.local_searcher.empty_ignored_nodes()
			else:
				self.local_searcher.fill_ignored_nodes(community)

		# # for debug
		# print(self.starting_nodes)

		if with_amend:
			self.amend_partition()

		return self.partition
