import networkx as nx
from LoadGraph import load_graph
from LocalSiwo import LocalSiwoCommunityDiscovery, LocalSiwoCommunityDetection
import time


def test_for_community_search(G):
	start = time.time()
	community_searcher = LocalSiwoCommunityDiscovery(G)
	for n in sorted(G.nodes()):
		community = community_searcher.community_search(start_node=n, with_amend=True)
		community_searcher.reset()
		print(n, ':\t', community, len(community))
	print(time.time() - start)


def test_for_community_detection(G):
	start = time.time()
	community_detector = LocalSiwoCommunityDetection(G)
	partition = community_detector.community_detection('random', overlap_enabled=False, with_amend=True)
	for e, p in enumerate(partition):
		print(e+1, p, len(p))
	print(time.time() - start)


if __name__ == '__main__':
	graph = load_graph('karate_edgelist.mtx')
	print(nx.info(graph))

	test_for_community_search(graph)
	print('-' * 50)

	test_for_community_detection(graph)
	print('-' * 50)
