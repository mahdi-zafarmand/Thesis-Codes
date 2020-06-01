file_address_source = 'dblp_edgelist1.txt'
file_address_destination = 'dblp_edgelist.txt'

nodes = dict()
list_of_edges = list()
counter = 0

with open(file_address_source, 'r') as file:
	lines = file.readlines()
	for line in lines:
		line = line.split()
		
		if not int(line[0]) in nodes:
			nodes[int(line[0])] = counter
			counter += 1

		if not int(line[1]) in nodes:
			nodes[int(line[1])] = counter
			counter += 1

		list_of_edges.append( (nodes[int(line[0])], nodes[int(line[1])]) )



with open(file_address_destination, 'w') as file:
	for edge in list_of_edges:
		to_write = str(edge[0]) + '\t' + str(edge[1]) + '\n'
		file.write(to_write)

print('Done')
