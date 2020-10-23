from CommunitySearch import CommunitySearcher
import networkx as nx
from copy import deepcopy


class ModularityLCommunityDiscovery(CommunitySearcher):
	def __init__(self, name):
		super().__init__(name)
