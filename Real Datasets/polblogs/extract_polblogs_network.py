import networkx as nx


gmlFile = 'polblogs.gml'

graph = nx.read_gml(gmlFile)
node2com = dict()
nodeMap1 = dict()
nodeMap2 = dict()
edgelist = list()

for e, node in enumerate(graph.nodes(data=True)):
	nodeMap1[node[0]] = e
	nodeMap2[e] = node[0]
	node2com[e] = node[1].get('value')

print('reading nodes = complete')

for n1, n2, edge in graph.edges(data=True):
	edgeToappend = (nodeMap1[n1], nodeMap1[n2])
	edgelist.append(edgeToappend)

print('reading edges = complete')

with open('polblogs_edgelist.mtx', 'w') as file:
	for edge in edgelist:
		lineTowrite = str(edge[0]) + '\t' + str(edge[1]) + '\n'
		file.write(lineTowrite)

with open('polblogs_groundtruth.mtx', 'w') as file:
	for node, com in node2com.items():
		lineTowrite = str(node) + '\t' + str(com) + '\n'
		file.write(lineTowrite)