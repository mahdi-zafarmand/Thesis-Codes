import fileinput

filename = 'karate_edgelist.mtx'

for line in fileinput.FileInput(filename, inplace=1):
    print(line.replace(' ', '\t'), end='')
