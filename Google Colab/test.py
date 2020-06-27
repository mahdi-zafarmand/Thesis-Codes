import local_siwo
import generate_lfr
import utils


network_specs = list()
network_specs.append({'n':100,      'tau1':20, 'tau2':10, 'mu':0.1, 'min_deg':10})
network_specs.append({'n':1000,     'tau1':20, 'tau2':10, 'mu':0.1, 'min_deg':10})
# network_specs.append({'n':10000,    'tau1':20, 'tau2':10, 'mu':0.1, 'min_deg':10})
# network_specs.append({'n':20000,    'tau1':20, 'tau2':10, 'mu':0.1, 'min_deg':10})

network_specs.append({'n':1000,    'tau1':20, 'tau2':10, 'mu':0.1, 'min_deg':10})
network_specs.append({'n':1000,    'tau1':20, 'tau2':10, 'mu':0.1, 'min_deg':20})
network_specs.append({'n':1000,    'tau1':20, 'tau2':10, 'mu':0.1, 'min_deg':40})
network_specs.append({'n':1000,    'tau1':20, 'tau2':10, 'mu':0.1, 'min_deg':80})

# for spec in network_specs:
#   network_path, groundtruth_path = generate_lfr.run_generator(num_nodes=spec['n'], tau1=spec['tau1'], tau2=spec['tau2'], mu=spec['mu'], min_degree=spec['min_deg'])
#   graph = utils.load_graph(network_path, weighted=False, self_loop=True)
#   local_siwo.community_detection(graph, groundtruth_path)
#   print('-' * 25)


network_path = 'Real Datasets/karate_edgelist.mtx'
groundtruth_path = 'Real Datasets/karate_groundtruth.mtx'
graph = utils.load_graph(network_path, weighted=False, self_loop=True)
local_siwo.community_search_for_all_nodes(graph, groundtruth_path)
