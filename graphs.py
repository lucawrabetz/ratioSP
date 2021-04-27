import random
import networkx as nx
import numpy as np

FILE_NUM = 4
DENSE_BASE_FILENAME = "_dense_.graph"
SPARSE_BASE_FILENAME = "_sparse_.graph"

# i is the size of the graph (n)
# dense graphs are erdos renyi with p = 0.8
# sparse graphs are erdos renyi with p = 0.5
stagger_val = 250
DENSE_FILE_NAMES = [str(i+stagger_val) +
                    DENSE_BASE_FILENAME for i in range(FILE_NUM)]
SPARSE_FILE_NAMES = [
    str(i+stagger_val) + SPARSE_BASE_FILENAME for i in range(FILE_NUM)]

# values for average and standard deviation for normal distribution for costs
mu_c = 250
mu_tau = 400
sigma = 50


def write_graph(filename, N, dense=True):
    connected = False
    if dense:
        prob = 0.8
    else:
        prob = 0.5

    while connected == False:
        node_node_matrix = [[0 for j in range(N)] for i in range(N)]
        exist_list = list(np.random.binomial(n=1, p=prob, size=N**2))
        cost_list = list(np.random.normal(mu_c, sigma, size=N**2))
        tau_list = list(np.random.normal(mu_tau, sigma, size=N**2))
        exists = 0
        counter = 0
        m = 0
        G = nx.Graph()

        for u in range(N):
            for v in range(u+1, N):
                exists = exist_list[counter]
                if exists:
                    G.add_edge(u, v, cost=int(
                        cost_list[counter]), tau=int(tau_list[counter]))
                    m += 1
                counter += 1

        # ensure graph is connected
        if nx.is_connected(G):
            connected = True
        else:
            continue

        # write edge list to file
        nx.write_edgelist(G, filename, data=['cost', 'tau'])

        # prepend n and m to first line
        with open(filename, 'r') as original:
            data = original.read()
        with open(filename, 'w') as new:
            new.write(str(N) + " " + str(m) + "\n" + data)

    return


for i in range(FILE_NUM):
    write_graph(DENSE_FILE_NAMES[i], i+stagger_val)
    write_graph(SPARSE_FILE_NAMES[i], i+stagger_val, False)
