import networkx as nx
import numpy as np
import random
import measures


num_nodes = 250             # number of nodes in the generated network
tau1 = 3.0                  # parameter for degree distribution: must be strictly greater than 1
tau2 = 1.5                  # parameter for community size distribution: must be strictly greater than 1
mu = 0.1                    # faction of links between communities: must be in [0, 0.5]
avg_degree = 10             # average degree of nodes: must be in [0, num_nodes]
max_degree = 20             # maximum degree of a node in the graph: must be less than num_ndoes
pdf = np.random.uniform     # the probability distribution function to generate probability values
prob_links_ratio = 0.1      # ratio of new probabilistic links which we add to the network
swap_ratio = 0.1            # ratio of number of low probability values that are swapped with high probability values to all probability values


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


# determine number of links of the network and number of uncertain links
num_links = G.number_of_edges()
num_uncertain_links = int(num_links * prob_links_ratio)
print('Number of all links in the network =', num_links)
print('Number of uncertain links with probabilty less than one =', num_uncertain_links)
print()


# create a list of probability values and sort them descendingly
probs = pdf(low=0.0, high=1.0, size=(num_uncertain_links, ))
probs = np.round(probs, decimals=2)
probs[::-1].sort()      
print('The distribution used to generate probabilities =', pdf.__name__)
print()


# determine the number of links inside communities (intra links) and between them (inter links)
num_intra_links = nx.algorithms.community.quality.intra_community_edges(G, comms)
num_inter_links = nx.algorithms.community.quality.inter_community_edges(G, comms)
print('Number of all links inside communities (intra links) =', num_intra_links)
print('Number of all links between communities (inter links) =', num_inter_links)
print()


# create probability arrays for intra links and inter links
num_swaps = int(swap_ratio * len(probs))
for i in random.sample(range(len(probs)), num_swaps):
    probs[i], probs[len(probs) - 1 - i] = probs[len(probs) - 1 - i], probs[i]

separator_index = int(num_intra_links / num_links * len(probs))
intra_probs = list(probs[:separator_index])
inter_probs = list(probs[separator_index:])


# updating the probability of edges
edges = list(G.edges(data=True))
random.shuffle(edges)
for edge in edges:
    n1, n2, info = edge
    if node2com[n1] == node2com[n2] and intra_probs:    # if edge is inside a community and there still exists a intra probability to assign
        p = intra_probs.pop()
        G.add_edge(n1, n2, prob=p)
    elif node2com[n1] != node2com[n2] and inter_probs:  # if edge is between communities and there still exists an inter probability to assign
        p = inter_probs.pop()
        G.add_edge(n1, n2, prob=p)
    else:                                               # if there is no proper probability value to assign (and in case of a certain link)
        G.add_edge(n1, n2, prob=1.0)


# create a list of edges and write them in the "edge_list.txt", 
# each edge is in a separate line.
edges = list(G.edges(data=True))
edges.sort()
with open('edge_list.txt', 'w') as file:
    for n1, n2, info in edges:
        prob = info.get('prob', 1.00)
        to_write = str(n1) + '\t' + str(n2) + '\t' + str(prob) + '\n'
        file.write(to_write)
    file.close()


# find the community index of a node and then write them in "node2com.txt", 
# the community index of each node is in a separate line.
with open('node2com.txt', 'w') as file:
    for node, com in node2com.items():
        to_write = str(com) + '\n'
        file.write(to_write)
    file.close()


# find the nodes of each community and write them in "com2nodes.txt",
# nodes of each community go to a separate line.
with open('com2nodes.txt', 'w') as file:
    for com, nodes in com2nodes.items():
        nodes.sort()
        to_write = ''
        for node in nodes:
            to_write += str(node) + ' '
        to_write = to_write[:-1] + '\n'
        file.write(to_write)
    file.close()


print('Done.')
