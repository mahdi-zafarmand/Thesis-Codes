import networkx as nx
import numpy as np
import random
import measures


num_nodes = 50              # number of nodes in the generated network
tau1 = 3.0                  # parameter for degree distribution: must be strictly greater than 1
tau2 = 1.5                  # parameter for community size distribution: must be strictly greater than 1
mu = 0.1                    # faction of links between communities: must be in [0, 0.5]
avg_degree = 10             # average degree of nodes: must be in [0, num_nodes]
max_degree = 20             # maximum degree of a node in the graph: must be less than num_ndoes


# create graph G via LFR, then determine its ground-truth communities
done = False
num_of_tries = 0
while not done and num_of_tries < 10:
    try:
        G = nx.LFR_benchmark_graph(num_nodes, tau1, tau2, mu, avg_degree, max_degree=max_degree)
        done = True
    except:
        num_of_tries += 1

if num_of_tries == 10:
    print("You need to change the parameters, this doesn't work.")
    exit(0)

comms = {frozenset(G.nodes[v]['community']) for v in G}
com2nodes = dict()
node2com = dict()
n = 0
for com in comms:
    com2nodes[n] = list(com)
    for node in com:
        node2com[node] = n
    n += 1

print(nx.info(G))
print('Modularity based on Newman formula =', round((measures.modularity(G, com2nodes)), 3))
print('Number of communities =', len(com2nodes))
for com_index, nodes in com2nodes.items():
    nodes.sort()
    print(com_index, '-> len of com =', len(nodes), '\t:', nodes)
print()

graph_name = 'synthetic_' + str(num_nodes) + '_'  + str(len(com2nodes)) + '.mtx'
nx.write_edgelist(G, graph_name, delimiter='\t', data=False)

groundTruth_name = graph_name[:-4] + '_ground_truth.mtx'
with open(groundTruth_name, 'w') as file:
	for node, com in sorted(node2com.items()):
		lineTowrite = str(node) + '\t' + str(com) + '\n'
		file.write(lineTowrite)