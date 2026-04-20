#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Dual graph visualization.

Public API
----------
plot_dual_graph(graph_id, output_path) -> bool
    Load graph_id from the DualEig/DualAdj library and save a PNG to output_path.
    Returns True on success, False if the graph_id is not found.
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")           # non-interactive backend — safe for web servers
import matplotlib.pyplot as plt
import networkx as nx

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)

from ClassesFunctions import loadEigenvalues, loadAdjMatrices


def _load_graph(graph_id: str):
    """Return the DualGraph object for graph_id, or None if not found."""
    try:
        n = int(graph_id.split("_")[0])
    except (ValueError, IndexError):
        return None

    eigen_file = os.path.join(_SCRIPT_DIR, "DualEig", f"{n}Eigen")
    adj_file   = os.path.join(_SCRIPT_DIR, "DualAdj",  f"V{n}AdjDG")

    if not os.path.isfile(eigen_file) or not os.path.isfile(adj_file):
        return None

    graphs = []
    loadEigenvalues(graphs, n, eigen_file)
    loadAdjMatrices(graphs, n, adj_file)

    for g in graphs:
        if g.graphID == graph_id:
            return g
    return None


def plot_dual_graph(graph_id: str, output_path: str) -> bool:
    """Generate a PNG visualization of graph_id and save it to output_path.

    Args:
        graph_id:    Dual graph ID string, e.g. '3_4' or '5_55'.
        output_path: Destination .png file path.

    Returns:
        True if the plot was saved successfully, False otherwise.
    """
    g = _load_graph(graph_id)
    if g is None:
        return False

    A = np.array(g.adjMatrix)
    n = len(A)

    # Build a MultiGraph so parallel edges and self-loops render correctly
    G = nx.MultiGraph()
    G.add_nodes_from(range(n))

    self_loop_nodes = set()
    for i in range(n):
        for j in range(i, n):
            count = int(A[i][j])
            if count == 0:
                continue
            if i == j:
                # self-loop (hairpin / loop helix)
                self_loop_nodes.add(i)
            else:
                for _ in range(count):
                    G.add_edge(i, j)

    # --- layout -------------------------------------------------------
    # Circular layout guarantees every node is equidistant, preventing
    # spring-layout artefacts where nearby nodes make parallel edges
    # between different pairs look like they share the same endpoints.
    if n == 1:
        pos = {0: np.array([0.0, 0.0])}
    else:
        pos = nx.circular_layout(G)

    # --- figure -------------------------------------------------------
    figsize = max(6, n * 1.4)          # larger canvas for bigger graphs
    fig, ax = plt.subplots(figsize=(figsize, figsize))
    ax.set_aspect("equal")
    ax.axis("off")

    # Edges between distinct nodes — draw each parallel edge with a different
    # curvature so they are visually separated (arrows=True enables connectionstyle).
    from collections import defaultdict
    edge_counts = defaultdict(int)
    for u, v in G.edges():
        key = (min(u, v), max(u, v))
        edge_counts[key] += 1

    drawn = defaultdict(int)
    for u, v in G.edges():
        key = (min(u, v), max(u, v))
        total = edge_counts[key]
        idx   = drawn[key]
        drawn[key] += 1
        if total == 1:
            rad = 0.0
        else:
            # Use a generous step so parallel edges are visually distinct
            # even when nodes sit close together in the spring layout.
            step = 0.45
            rad  = -step * (total - 1) / 2 + idx * step
        nx.draw_networkx_edges(
            G, pos, ax=ax,
            edgelist=[(u, v)],
            edge_color="#c0392b",
            width=3.0,
            alpha=0.9,
            arrows=True,
            arrowstyle="-",
            connectionstyle=f"arc3,rad={rad:.3f}",
            min_source_margin=20,
            min_target_margin=20,
        )

    # Self-loops: draw as circles near each self-loop node
    loop_radius = 0.14
    loop_offset = 0.22
    for node in self_loop_nodes:
        x, y = pos[node]
        loop = plt.Circle((x, y + loop_offset), loop_radius,
                           color="#c0392b", fill=False, linewidth=3.5)
        ax.add_patch(loop)

    # Expand axis limits so nodes and self-loop patches are never clipped
    ax.autoscale_view()
    xl, yl = ax.get_xlim(), ax.get_ylim()
    pad_xy = 0.35
    top_extra = (loop_offset + loop_radius + 0.05) if self_loop_nodes else pad_xy
    ax.set_xlim(xl[0] - pad_xy, xl[1] + pad_xy)
    ax.set_ylim(yl[0] - pad_xy, yl[1] + top_extra)

    # Nodes
    nx.draw_networkx_nodes(
        G, pos, ax=ax,
        node_color="#c0392b",
        node_size=900,
        edgecolors="#922b21",
        linewidths=2,
    )

    # Labels  (1-indexed to match the dual graph convention)
    labels = {i: str(i + 1) for i in range(n)}
    nx.draw_networkx_labels(
        G, pos, labels, ax=ax,
        font_size=15,
        font_color="white",
        font_weight="bold",
    )

    ax.set_title(
        f"Dual Graph  {graph_id}\n{n} {'vertex' if n == 1 else 'vertices'}",
        fontsize=16, fontweight="bold", pad=16,
    )

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return True


if __name__ == "__main__":
    # CLI usage: python plot_dual_graph.py 3_4 output.png
    import argparse
    parser = argparse.ArgumentParser(description="Plot a dual graph topology.")
    parser.add_argument("graph_id", help="Dual graph ID, e.g. 3_4")
    parser.add_argument("output",   help="Output PNG path")
    args = parser.parse_args()

    ok = plot_dual_graph(args.graph_id, args.output)
    if ok:
        print(f"Saved: {args.output}")
    else:
        print(f"Graph '{args.graph_id}' not found in library.")
        sys.exit(1)
