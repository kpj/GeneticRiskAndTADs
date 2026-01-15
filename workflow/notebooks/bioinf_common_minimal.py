"""Some functions taken from https://github.com/kpj/bioinf_common/."""

import sys
import random
from typing import Any
from collections.abc import Sequence, Mapping

import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests

import networkx as nx

import seaborn as sns
import matplotlib.pyplot as plt


def multipletests_nan(
    pval_list: Sequence[float], method: str = "fdr_bh"
) -> Sequence[float]:
    """Multiple testing correction on lists with NaNs."""
    pval_list = np.asarray(pval_list)

    (nan_idx,) = np.where(np.isnan(pval_list))
    pval_list_nonan = pval_list[~np.isnan(pval_list)]

    if len(pval_list_nonan) == 0:
        assert np.isnan(pval_list).all()  # contains only np.nan
        pval_corr = pval_list
    else:
        _, pval_corr, _, _ = multipletests(pval_list_nonan, method=method)
        for i in nan_idx:
            pval_corr = np.insert(pval_corr, i, np.nan)

    return pval_corr


def add_identity(ax, *line_args, **line_kwargs):
    """Add identity (y=x) line to given axes."""
    # line plot
    (identity,) = ax.plot([], [], *line_args, **line_kwargs)

    # react to layout/data changes
    def callback(ax):
        low_x, high_x = ax.get_xlim()
        low_y, high_y = ax.get_ylim()
        low = max(low_x, low_y)
        high = min(high_x, high_y)
        identity.set_data([low, high], [low, high])

    callback(ax)
    ax.callbacks.connect("xlim_changed", callback)
    ax.callbacks.connect("ylim_changed", callback)

    return ax


def get_fraction_of_nonisolated_nodes(
    nodes: Sequence[Any], graph: nx.Graph, nc_type: str = "overall"
) -> float:
    """Control as
    * len(non_isolated_nodes) / len(all_nodes)
    * len(non_isolated_nodes) / (len(isolated_nodes) + 1)
    """
    ns = set(nodes)
    isolated_nodes = []
    non_isolated_nodes = []
    for n in nodes:
        neighs = graph.neighbors(n)
        if (ns - {n}).intersection(neighs):
            non_isolated_nodes.append(n)
        else:
            isolated_nodes.append(n)

    assert len(isolated_nodes) == len(set(isolated_nodes))
    assert len(non_isolated_nodes) == len(set(non_isolated_nodes))
    assert len(isolated_nodes) <= len(nodes)
    assert len(non_isolated_nodes) <= len(nodes)

    if nc_type == "overall":
        return len(non_isolated_nodes) / len(nodes)
    elif nc_type == "overiso":
        return len(non_isolated_nodes) / (len(isolated_nodes) + 1)
    else:
        raise RuntimeError(f'Invalid NC-type "{nc_type}"')


def compute_network_coherence(
    graph: nx.Graph,
    nodes: Sequence[Any],
    random_nodes: Sequence[Any] | None = None,
    reps: int = 1000,
    fail_on_zero_std: bool = False,
    nc_type: str = "overall",
) -> float:
    """Compute z-score of non-isolated node fraction (network coherence)"""
    if len(nodes) == 0:
        return np.nan

    assert len(set(nodes)) == len(nodes), nodes
    random_nodes = list(random_nodes or graph.nodes())

    rand_node_list: Sequence[Any] = []
    while len(rand_node_list) < reps:
        cur = random.sample(random_nodes, len(nodes))

        assert len(set(cur)) == len(nodes)
        rand_node_list.append(cur)

    # network coherence
    act_num = get_fraction_of_nonisolated_nodes(nodes, graph, nc_type)
    rand_nums = [
        get_fraction_of_nonisolated_nodes(rn, graph, nc_type) for rn in rand_node_list
    ]

    # handle case where no z-score can be computed
    if np.std(rand_nums) == 0:
        if fail_on_zero_std or len(nodes) == 1:
            return np.nan
        else:  # retry if sensible
            next_rep_num = reps * 3
            print(
                "Uniform random results, retrying with more iterations "
                f"({reps} -> {next_rep_num})",
                file=sys.stderr,
            )
            return compute_network_coherence(
                graph,
                nodes,
                random_nodes=random_nodes,
                reps=next_rep_num,
                fail_on_zero_std=True,
            )

    z_score = (act_num - np.mean(rand_nums)) / np.std(rand_nums)
    return z_score


def annotated_barplot(anno_kws: Mapping[str, Any] | None = None, **kwargs: Any) -> None:
    anno_kws = anno_kws or {}

    label_offset = anno_kws.pop("label_offset", 20)
    label_size = anno_kws.pop("label_size", plt.rcParams["font.size"])
    if len(anno_kws) > 0:
        raise TypeError(f"Invalid parameters: {anno_kws.keys()}")

    g = sns.barplot(**kwargs)

    # for i, row in enumerate(kwargs['data'].itertuples()):
    #     g.annotate(
    #         f'{row.value:,.1f}', (i, row.value), xycoords='data',
    #         ha='center', xytext=(0, 3), textcoords='offset pixels')

    for p in g.patches:
        # handle custom baseline
        height = p.get_height() + kwargs.get("bottom", 0)

        g.annotate(
            f"{height:,.1f}",
            (p.get_x() + p.get_width() / 2.0, height),
            ha="center",
            va="center",
            color="gray",
            fontsize=label_size,
            xytext=(0, label_offset),
            textcoords="offset pixels",
        )
