import networkx as nx
import numpy as np
import random
import measures
import time


num_nodes = 20000            	# number of nodes in the generated network
tau1 = 20                  		# parameter for degree distribution: must be strictly greater than 1
tau2 = 10                  		# parameter for community size distribution: must be strictly greater than 1
mu = 0.1                    	# faction of links between communities: must be in [0, 0.5]
min_degree = 50 		        
max_degree = num_nodes
num_of_tries = 10				# number of tries to generate the network



def generate_the_network(num_nodes, tau1, tau2, mu, min_degree, max_degree):
	i = 0
	done = False
	while i < num_of_tries and done == False:
		print('i =', i)
		try:
		    G = nx.LFR_benchmark_graph(num_nodes, tau1, tau2, mu, min_degree=min_degree, max_degree=max_degree)
		    done = True
		except Exception as e:
			pass
		i += 1

	if done == False:
		return None

	return G


def extract_comms(G):
	comms = {frozenset(G.nodes[v]['community']) for v in G}
	com2nodes = dict()
	node2com = dict()
	n = 0
	for com in comms:
	    com2nodes[n] = list(com)
	    for node in com:
	        node2com[node] = n
	    n += 1
	return node2com, com2nodes


def print_graph_info(G, com2nodes):
	print(nx.info(G))
	print('Modularity based on Newman formula =', round((measures.modularity(G, com2nodes)), 3))
	print('Number of communities =', len(com2nodes))

	print('sizes of comms = [', end="")
	for i in range(len(com2nodes)):
		if i == len(com2nodes) - 1:
			print(len(com2nodes[i]), end="")
		else:
			print(len(com2nodes[i]), ', ', end="")
	print(']\n')


def write_to_file(G, node2com, com2nodes):
	num_nodes = G.number_of_nodes()
	graph_name = 'synthetic_' + str(num_nodes) + '_'  + str(len(com2nodes)) + '.mtx'
	nx.write_edgelist(G, graph_name, delimiter='\t', data=False)

	groundTruth_name = graph_name[:-4] + '_ground_truth.mtx'
	with open(groundTruth_name, 'w') as file:
		for node, com in sorted(node2com.items()):
			lineTowrite = str(node) + '\t' + str(com) + '\n'
			file.write(lineTowrite)


def _main():
	start_time = time.time()
	synthetic_network = generate_the_network(num_nodes, tau1, tau2, mu, min_degree, max_degree)
	if synthetic_network != None:
		node2com, com2nodes = extract_comms(synthetic_network)

		print_graph_info(synthetic_network, com2nodes)
		write_to_file(synthetic_network, node2com, com2nodes)

	finish_time = time.time()
	print('\nDone in %.4f seconds.' %(finish_time - start_time))

		

if __name__ == '__main__':
	_main()