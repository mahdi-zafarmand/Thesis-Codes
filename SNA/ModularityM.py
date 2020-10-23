from CommunitySearch import CommunitySearcher
import networkx as nx
from copy import deepcopy


class ModularityMCommunityDiscovery(CommunitySearcher):
	def __init__(self, name):
		super().__init__(name)
