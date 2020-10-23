from copy import deepcopy
import utils


class CommunityDetector:
	# minimum_improvement = 0.000001    # I think it is not practical.

	def __init__(self, name, graph):
		# initialize the object.
		self.name = name
		self.graph = deepcopy(graph)
		self.local_searcher = None
		self.nodes_to_be_ignored = []
		self.overlap_enable = False
		self.node_selection_mode = 'random'
		self.starting_nodes = []
		self.discovered_nodes = set()
		self.number_discovered_nodes = 0
		self.partition = []
		self.remove_self_loops()

	def remove_self_loops(self):
		# algorithms tend to work better if there is no self-loop in the given graph, so we call this method at first.
		utils.remove_self_loops(self.graph)

	def update_discovered_nodes_info(self):
		# after a community is found and appended to the 'self.partition', we need to update 2 attributes of the class.
		for node in self.partition[-1]:
			if (node in self.discovered_nodes) is False:
				self.discovered_nodes.add(node)
				self.number_discovered_nodes += 1

	def determine_starting_node(self, mode):
		# check the validity of the given mode, sets the 'self.node_selection_mode'
		# and adds the newest starting node to 'self.starting_nodes'.
		if mode == 'random':
			self.node_selection_mode = 'random'
			self.starting_nodes.append(self.find_random_node_in_graph())
		elif mode == 'highest_degree':
			self.node_selection_mode = 'highest_degree'
			self.starting_nodes.append(self.find_highest_degree_node_in_graph())
		else:
			print("Invalid mode to select the starting node!", end=' ')
			print("It should be either 'random' or 'highest_degree'.")
			exit(-1)

	def find_random_node_in_graph(self):
		# returns a random node among nodes of the 'self.graph', except the ones that should be ignored.
		all_nodes_to_be_ignore = [self.nodes_to_be_ignored, self.discovered_nodes]
		return utils.find_random_node_in_graph(self.graph, all_nodes_to_be_ignore)

	def find_highest_degree_node_in_graph(self):
		# returns the node with the highest degree among nodes of the 'self.graph', except the ones that should be ignored.
		all_nodes_to_be_ignore = [self.nodes_to_be_ignored, self.discovered_nodes]
		return utils.find_highest_degree_node_in_graph(self.graph, all_nodes_to_be_ignore)

	def community_detection(self, node_selection_mode, overlap_enabled=False, with_amend=True):
		# THE MAIN FUNCTION OF THE CLASS, finds all communities of the graph.
		while self.number_discovered_nodes < self.graph.number_of_nodes():
			self.determine_starting_node(node_selection_mode)
			community = self.local_searcher.community_search(start_node=self.starting_nodes[-1], with_amend=with_amend)
			self.partition.append(community)
			self.local_searcher.reset()
			self.update_discovered_nodes_info()

			if overlap_enabled:
				self.local_searcher.empty_ignored_nodes()
			else:
				self.local_searcher.fill_ignored_nodes(community)

		if with_amend:
			self.amend_partition()
		return self.partition

	def compute_modularity(self, neighbor_node):
		# may not be required.
		pass

	def amend_partition(self):
		# as most methods works based on triangles, we should remedy the results for communities smaller than size 3.
		communities_of_size_one = [community for community in self.partition if len(community) == 1]
		communities_of_size_two = [community for community in self.partition if len(community) == 2]

		for community in communities_of_size_one:
			self.partition.remove(community)
		for community in communities_of_size_two:
			self.partition.remove(community)

		if len(communities_of_size_one) > 0:
			self.amend_partition_for_size_one(communities_of_size_one)
		if len(communities_of_size_two) > 0:
			self.amend_partition_for_size_two(communities_of_size_two)

	def amend_partition_for_size_one(self, communities_of_size_one):
		# deals with the communities containing only one node, if 'with_amend' = True.
		for community in communities_of_size_one:
			neighbors = set(self.graph.neighbors(community[0]))
			strength_dict = {}
			for neighbor in neighbors:
				for i in range(len(self.partition)):
					if neighbor in self.partition[i] and self.graph.has_edge(community[0], neighbor):
						# next line is different than what is was in previous code: 'weight' -> 'strength'
						strength_dict[i] = strength_dict.get(i, 0.0) + self.graph[community[0]][neighbor].get('strength', 0.0)
						if not self.overlap_enable:
							break
			self.amend_partition_helper(community, strength_dict)

	def amend_partition_for_size_two(self, communities_of_size_two):
		# deals with the communities containing only two nodes, if 'with_amend' = True.
		for community in communities_of_size_two:
			neighbors = set(self.graph.neighbors(community[0]))
			neighbors.update(self.graph.neighbors(community[1]))
			strength_dict = {}
			for neighbor in neighbors:
				for i in range(len(self.partition)):
					if neighbor in self.partition[i]:
						if self.graph.has_edge(community[0], neighbor):
							# next line is different than previous code: 'weight' -> 'strength' : seems OK
							strength_dict[i] = strength_dict.get(i, 0.0) + self.graph[community[0]][neighbor].get('strength', 0.0)
						if self.graph.has_edge(community[1], neighbor):
							# next line is different than previous code: 'weight' -> 'strength' : seems OK
							strength_dict[i] = strength_dict.get(i, 0.0) + self.graph[community[1]][neighbor].get('strength', 0.0)
						if not self.overlap_enable:
							break
			self.amend_partition_helper(community, strength_dict)

	def amend_partition_helper(self, community, strength_dict):
		# a helper method, amends the partition by merging communities of sizes 1 and 2 into other proper communities.
		index_best_community_to_merge_into = list(strength_dict.keys())[0]
		for index_community in strength_dict:
			if strength_dict[index_community] > strength_dict[index_best_community_to_merge_into]:
				index_best_community_to_merge_into = index_community
		self.partition[index_best_community_to_merge_into].extend(community)

		# next line to sort the community is only for a better representation, can be ignored to boost performance.
		self.partition[index_best_community_to_merge_into].sort()
