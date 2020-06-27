import networkx as nx


gmlFile = 'polbooks.gml'

graph = nx.read_gml(gmlFile)
node2com = dict()
nodeMap1 = dict()
nodeMap2 = dict()
edgelist = list()

for e, node in enumerate(graph.nodes(data=True)):
	nodeMap1[node[0]] = e
	nodeMap2[e] = node[0]
	if node[1].get('value') == 'n':
		node2com[e] = 0
	elif node[1].get('value') == 'c':
		node2com[e] = 1
	else:
		node2com[e] = 2


for n1, n2, edge in graph.edges(data=True):
	edgeToappend = (nodeMap1[n1], nodeMap1[n2])
	edgelist.append(edgeToappend)

with open('polbooks_edgelist.mtx', 'w') as file:
	for edge in edgelist:
		lineTowrite = str(edge[0]) + '\t' + str(edge[1]) + '\n'
		file.write(lineTowrite)

with open('polbooks_groundtruth.mtx', 'w') as file:
	for node, com in node2com.items():
		lineTowrite = str(node) + '\t' + str(com) + '\n'
		file.write(lineTowrite)