import networkx as nx
import utils
import time
from copy import deepcopy


_N_STEP = 4
_N_ITER = 3


def insise_triangles(graph, community):
    count = 0
    for node1 in community:
        for node2 in community:
            if node1 != node2:
                for node3 in community:
                    if node3 != node1 and node3 != node2:
                        if graph.has_edge(node1, node2) and graph.has_edge(node1, node3) and graph.has_edge(node2, node3):
                            count += 1
    return int(count / 6)


def outside_triangles(graph, community):
    count = 0
    for node in community:
        for neigh1 in graph.neighbors(node):
            for neigh2 in graph.neighbors(node):
                if neigh1 != neigh2 and not neigh1 in community and not neigh2 in community:
                    if graph.has_edge(neigh1, neigh2):
                        count += 1
    return int(count / 2)


def calc_T_in(graph, node, community, Tin_prev):
    increment = 0.0
    neighbors = list(graph.neighbors(node))

    for neigh1 in neighbors:
        for neigh2 in neighbors:
            if neigh1 != neigh2 and graph.has_edge(neigh1, neigh2):
                if neigh1 in community and neigh2 in community:
                    increment += 1

    return Tin_prev + 0.5 * increment


def calc_T_out(graph, node, community, Tout_prev):
    increment = 0.0
    decrement = 0.0
    neighbors = list(graph.neighbors(node))

    for neigh1 in neighbors:
        for neigh2 in neighbors:
            if neigh1 != neigh2 and graph.has_edge(neigh1, neigh2):
                if not neigh1 in community and not neigh2 in community:
                    increment += 1
                elif neigh1 in community and not neigh2 in community:
                    decrement += 1

    return Tout_prev + 0.5 * increment - decrement


def calc_T(Tin, Tout):
    if Tin > Tout:
        return Tin * (Tin - Tout)
    return 0.0


def find_next_neighbor(graph, node, depth=2):
    if depth < _N_STEP:
        the_neighbors = dict()
        for neigh in graph.neighbors(node):
            the_neighbors[neigh] = find_next_neighbor(graph, neigh, depth + 1)

    else:
        the_neighbors = list(graph.neighbors(node))
        the_neighbors.sort()

    return the_neighbors


def read_trajectory(trajectory):
    type_values = type(list(trajectory.values())[0])
    if type_values == dict:
        output = dict()
        for k,v in trajectory.items():
            output[k] = read_trajectory(trajectory[k])
        return output
    elif type_values == list:
        output = list()
        for k, v in trajectory.items():
            for v_i in v:
                output.append((k, v_i))
        return output
    else:
        print('ERROR, THIS SHOULD NOT HAPPEN!')
        exit(0)


def clear_invalid_states(forward_states, community):
    i = 0
    while i < len(forward_states):
        if len(set(forward_states[i])) != len(forward_states[i]):
            forward_states.pop(i)
            i -= 1
        elif len(set(forward_states[i]).intersection(community)) > 0:
            forward_states.pop(i)
            i -= 1          
        i += 1
    return forward_states


def find_forward_states(graph, community):
    trajectory = dict()
    neighbors = list()
    for node in community:
        neighbors.extend(graph.neighbors(node))
    neighbors = list(set(neighbors))

    for neigh in sorted(neighbors):
        trajectory[neigh] = find_next_neighbor(graph, neigh)

    states = read_trajectory(trajectory)
    while type(states) != list:
        states_copy = deepcopy(states)
        states = read_trajectory(states_copy)

    forward_states = list()
    temp = list()
    for state in states:
        while type(state[1]) != int:
            temp.append(state[0])
            state = state[1]
        temp.append(state[0])
        temp.append(state[1])
        forward_states.append(temp[:])
        temp.clear()

    return clear_invalid_states(forward_states, community)


def compute_forward_scores(graph, Tin, Tout, community, forward_states):
    forward_scores = dict()
    inward_scores = dict()
    outward_scores = dict()

    for state in forward_states:
        Tin_state = Tin
        Tout_state = Tout
        for next_node in state:
            Tin_state = calc_T_in(graph, next_node, community, Tin_state)
            Tout_state = calc_T_out(graph, next_node, community, Tout_state)
            community.append(next_node)

        for i in range(_N_STEP):
            community.pop()

        inward_scores[tuple(state)] = Tin_state
        outward_scores[tuple(state)] = Tout_state
        forward_scores[tuple(state)] = calc_T(Tin_state, Tout_state)

    return forward_scores, inward_scores, outward_scores


def find_best_state(forward_scores, external_scores):
    states = list(forward_scores.keys())

    best_state = states[0]
    for i in range(1, len(states)):
        if forward_scores[states[i]] > forward_scores[states[i]]:
            best_state = states[i]
        elif forward_scores[states[i]] == forward_scores[states[i]]:
            if external_scores[states[i]] < external_scores[best_state]:
                best_state = states[i]

    return best_state


def find_best_next_node(graph, community, shell_set, Tin_prev, Tout_prev):
    dict_of_Tin = dict()
    dict_of_Tex = dict()
    dict_of_T = dict()

    for node in shell_set:
        T_in_node = calc_T_in(graph, node, community, Tin_prev)
        T_out_node = calc_T_out(graph, node, community, Tout_prev)
        dict_of_Tin[node] = T_in_node
        dict_of_Tex[node] = T_out_node
        dict_of_T[node] = calc_T(T_in_node, T_out_node)

    best_node = shell_set[0]
    for i in range(1, len(shell_set)):
        if dict_of_T[shell_set[i]] > dict_of_T[best_node]:
            best_node = shell_set[i]
        elif dict_of_T[shell_set[i]] == dict_of_T[best_node]:
            if dict_of_Tex[shell_set[i]] < dict_of_Tex[best_node]:
                best_node = shell_set[i]

    return best_node, dict_of_Tin[best_node], dict_of_Tex[best_node]



def community_search(graph, intended_node):
    community = list()
    shell_set = list()

    community.append(intended_node)
    Tin = calc_T_in(graph, intended_node, [], 0.0)
    Tout = calc_T_out(graph, intended_node, [], 0.0)

    num_iter = 0

    while len(community) < graph.number_of_nodes() and num_iter != _N_ITER:
        num_iter += 1
        old_score = calc_T(Tin, Tout)

        forward_states = find_forward_states(graph, community)
        forward_scores, internal_scores, external_scores = compute_forward_scores(graph, Tin, Tout, community, forward_states)
        best_next_node = find_best_state(forward_scores, external_scores)[0]

        Tin_new = calc_T_in(graph, best_next_node, community, Tin)
        Tout_new = calc_T_out(graph, best_next_node, community, Tout)
        new_score = calc_T(Tin_new, Tout_new)

        if new_score >= old_score:
            Tin, Tout = Tin_new, Tout_new
            old_score = new_score
            community.append(best_next_node)
            # print('Tin =', Tin, ' -> ', insise_triangles(graph, community), '\t\t', 'Tout =', Tout, ' -> ', outside_triangles(graph, community), '\n')    # uncomment only to check
        else:
            break
    print('END OF BEAM SEARCH PHASE!')

    shell_set = list()
    for node in community:
        shell_set.extend(graph.neighbors(node))
    shell_set = list(set(shell_set))

    while len(community) < graph.number_of_nodes():
        old_score = calc_T(Tin, Tout)
        best_node, Tin_new, Tout_new = find_best_next_node(graph, community, shell_set, Tin, Tout)
        new_score = calc_T(Tin_new, Tout_new)
        if new_score >= old_score:
            Tin, Tout = Tin_new, Tout_new
            new_neighbors = list(set(graph.neighbors(best_node)) - set(community))    
            community.append(best_node)
            shell_set.extend(new_neighbors)
            shell_set = list(set(shell_set))
            shell_set.remove(best_node)

            if not shell_set:
                break
        
        else:
            break


    return community


def _main():
    start_time = time.time()

    args = utils.create_argument_parser()
    graph = utils.load_graph(args.dataset, False)

    intended_node = int(args.output)

    community = community_search(graph, intended_node)
    print('community =', community, len(community))

    finish_time = time.time()
    print('\nDone in %.4f seconds.' %(finish_time - start_time))


if __name__ == "__main__":
    _main()