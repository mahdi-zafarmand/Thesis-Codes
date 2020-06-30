import generate_lfr
import utils
import local_siwo
import local_modularity_r
import local_modularity_m
import local_modularity_l
import local_triads
import local_nstep
import os


network_specs = list()
# # 1. if the number of nodes increases (100, 500, 1K, 5K, 10K), run all methods and compare (considering Q, NMI, and run time).
# network_specs.append({'n':100,      'tau1':20, 'tau2':10, 'mu':0.1, 'min_deg':10})
# network_specs.append({'n':500,      'tau1':20, 'tau2':10, 'mu':0.1, 'min_deg':10})
# network_specs.append({'n':1000,     'tau1':20, 'tau2':10, 'mu':0.1, 'min_deg':10})
# network_specs.append({'n':5000,     'tau1':20, 'tau2':10, 'mu':0.1, 'min_deg':10})
# network_specs.append({'n':10000,    'tau1':20, 'tau2':10, 'mu':0.1, 'min_deg':10})

# # 2. if the density of nodes increases (for a network with |V| = 1000), run all methods and compare (considering Q, NMI, and run time).
# network_specs.append({'n':1000,    'tau1':20, 'tau2':10, 'mu':0.1, 'min_deg':10})
# network_specs.append({'n':1000,    'tau1':20, 'tau2':10, 'mu':0.1, 'min_deg':20})
# network_specs.append({'n':1000,    'tau1':20, 'tau2':10, 'mu':0.1, 'min_deg':40})
# network_specs.append({'n':1000,    'tau1':20, 'tau2':10, 'mu':0.1, 'min_deg':60})
# network_specs.append({'n':1000,    'tau1':20, 'tau2':10, 'mu':0.1, 'min_deg':80})
# network_specs.append({'n':1000,    'tau1':20, 'tau2':10, 'mu':0.1, 'min_deg':100})

# # 3. if the mixing parameter decreases (for a network with |V| = 1000), run all methods and compare (considering Q, NMI, and run time).
# network_specs.append({'n':1000,     'tau1':20, 'tau2':10, 'mu':0.10, 'min_deg':10})
# network_specs.append({'n':1000,     'tau1':20, 'tau2':10, 'mu':0.15, 'min_deg':10})
# network_specs.append({'n':1000,     'tau1':20, 'tau2':10, 'mu':0.20, 'min_deg':10})
# network_specs.append({'n':1000,     'tau1':20, 'tau2':10, 'mu':0.25, 'min_deg':10})
# network_specs.append({'n':1000,     'tau1':20, 'tau2':10, 'mu':0.30, 'min_deg':10})
# network_specs.append({'n':1000,     'tau1':20, 'tau2':10, 'mu':0.35, 'min_deg':10})
# network_specs.append({'n':1000,     'tau1':20, 'tau2':10, 'mu':0.40, 'min_deg':10})


# # generate networks
# for i in range(len(network_specs)):
# 	spec = network_specs[i]
# 	network_path, groundtruth_path = generate_lfr.run_generator(num_nodes=spec['n'], tau1=spec['tau1'], tau2=spec['tau2'], mu=spec['mu'], min_degree=spec['min_deg'], prefix=str(i+1))
# 	# graph = utils.load_graph(network_path, weighted=False, self_loop=True)
# 	# print('Number of detected communities =', len(local_siwo.community_detection(graph, groundtruth_path)))
# 	print('-' * 25)


# some small tests:
# network_path = 'Synthetic Datasets/synthetic_250_20.mtx'
# groundtruth_path = 'Synthetic Datasets/synthetic_250_20_ground_truth.mtx'
# print(network_path,'\n')
# graph = utils.load_graph(network_path, weighted=False, self_loop=True)
# local_triads.community_search_for_all_nodes(graph, groundtruth_path)
# print('Number of detected communities =', len(local_triads.community_detection(graph, groundtruth_path)))


# read networks
directory_in_str = 'experiment_'
experiments = dict()
for i in range(1, 4):
	experiments[i] = list()
	directory = os.fsencode(directory_in_str + str(i))
	temp_dict = dict()
	for file in sorted(os.listdir(directory)):
		filename = os.fsdecode(file)
		if filename[2] == 'c':
			temp_dict['c'] = filename
		elif filename[2] == 'e':
			temp_dict['e'] = filename
			experiments[i].append(dict(temp_dict))
			temp_dict.clear()

# methods = [local_modularity_r, local_modularity_m, local_modularity_l, local_triads, local_nstep, local_siwo]
methods = [local_modularity_r, local_modularity_m, local_modularity_l, local_triads, local_siwo]
for exp_num, exp_args in experiments.items():
	for exp_ins in exp_args:
		network_path = 'experiment_' + str(exp_num) + '/' + exp_ins['e']
		groundtruth_path = 'experiment_' + str(exp_num) + '/' + exp_ins['c']
		print(network_path, ':')
		for method in methods:
			graph = utils.load_graph(network_path, weighted=False, self_loop=True)
			try:
				method.community_search_for_all_nodes(graph, groundtruth_path)
			except expression as identifier:
				print('SOME PROBLEM WITH METHOD', str(method), 'FOR COMMUNITY SEARCH')
			
			try:
				print('Number of detected communities using method="', str(method), '" =', len(method.community_detection(graph, groundtruth_path)))
			except expression as identifier:
				print('SOME PROBLEM WITH METHOD', str(method), 'FOR COMMUNITY DETECTION')

			print('----------------------------------------------------------------------------------------------------------------------------')
		print()
	print('**************************************************\n**************************************************')