from __future__ import print_function

import networkx as nx


def checkIsland(node, link, casename):
    link_edges = []
    # Create a list of all the edges
    for link_ele in link:
        link_edges.append((link_ele.from_bus, link_ele.to_bus))
    graph = nx.from_edgelist(link_edges)
    islands = list(nx.connected_components(graph))
    print('Total number of islands in this system is %d' % len(islands))
    return islands, graph


def linkNodePhaseCompatibility(link, node):
    pass
