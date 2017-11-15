"""
Store more elaborate measure computations
"""

import sys
import random

import numpy as np


def get_fraction_of_nonisolated_nodes(nodes, graph):
    """ Read function name
    """
    ns = set(nodes)
    non_isolated_nodes = []
    for n in nodes:
        neighs = graph.neighbors(n)
        if (ns-{n}).intersection(neighs):
            non_isolated_nodes.append(n)

    assert len(non_isolated_nodes) == len(set(non_isolated_nodes))
    assert len(non_isolated_nodes) <= len(nodes)

    return len(non_isolated_nodes) / len(nodes)

def compute_network_coherence(graph, nodes, random_nodes=None, reps=1000, fail_on_zero_std=False):
    """ Compute z-score of non-isolated node fraction (network coherence)
    """
    assert len(set(nodes)) == len(nodes)
    random_nodes = random_nodes or graph.nodes()

    rand_node_list = []
    while len(rand_node_list) < reps:
        cur = random.sample(random_nodes, len(nodes))

        assert len(set(cur)) == len(nodes)
        rand_node_list.append(cur)

    # network coherence
    g = lambda ns: get_fraction_of_nonisolated_nodes(ns, graph)
    act_num = g(nodes)
    rand_nums = [g(rn) for rn in rand_node_list]

    # handle case where no z-score can be computed
    if np.std(rand_nums) == 0:
        if fail_on_zero_std or len(nodes) == 1:
            return np.nan
        else:  # retry if sensible
            next_rep_num = reps*3
            print(
                'Uniform random results, retrying with more iterations '
                f'({reps} -> {next_rep_num})', file=sys.stderr)
            return compute_network_coherence(
                graph, nodes, random_nodes=random_nodes,
                reps=next_rep_num, fail_on_zero_std=True)

    z_score = (act_num - np.mean(rand_nums)) / np.std(rand_nums)
    return z_score
