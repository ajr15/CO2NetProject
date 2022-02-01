# script to analyze different variation of minimal vertex cover algorithms and analyses them
import os
import networkx as nx
from networkx.algorithms.approximation import vertex_cover
import numpy as np
import random
from ..commons import get_graph, distance_distribution
import settings


def project_to_type(G, target_type):
    new_G = nx.Graph()
    for node in G:
        if type(node) is target_type:
            neighbors = [nn for n in nx.all_neighbors(G, node) for nn in nx.all_neighbors(G, n)]
            for n in neighbors:
                new_G.add_edge(node, n)
    return new_G

def networkx_mvc_algorithm(G, weight_func=None):
    random.seed(0)
    np.random.seed(0)
    nG = project_to_type(G, str)
    if not weight_func is None:
        nx.classes.function.set_node_attributes(nG, {n: weight_func(n) for n in nG}, "cost")
        ws = "cost"
    else:
        ws = None
    return vertex_cover.min_weighted_vertex_cover(nG, ws)


def total_degree(rxn_graph, G, specie):
    """Greedy metric that calculates total degree of specie"""
    return len([s for s in nx.all_neighbors(G, specie.identifier)])

def percolation_degree(rxn_graph, G, specie):
    before = len(rxn_graph.reactions)
    after = len(rxn_graph.remove_specie(specie).reactions)
    return before - after

def greedy_mvc(rxn_graph, greedy_metric: callable):
    """Greedy algorithm for finding the MVC of a reaction graph.
    ARGS:
        - rxn_graph (TorinaNet.RxnGraph): reaction graph for analysis
        - greedy_metric (callable): function to calculate the score for choosing node in the greedy algorithm
                                    greedy metric is called as f(rxn_graph, networkx_rxn_graph, specie)"""
    # init
    _rxn_graph = rxn_graph.copy()
    mvc = set()
    while True:
        # randomly select reaction (all are uncovered)
        rxn = np.random.choice(_rxn_graph.reactions)
        # add to mvc product / reactant with greedy metric value
        G = _rxn_graph.to_networkx_graph()
        max_x = 0
        specie = ""
        for s in rxn.reactants + rxn.products:
            x = greedy_metric(_rxn_graph, G, s)
            if x > max_x and not s in _rxn_graph.reactant_species:
                max_x = x
                specie = s
        mvc.add(specie._get_id_str())
        print(len(_rxn_graph.reactions))
        # clear the network from dependent reactions
        _rxn_graph = _rxn_graph.remove_specie(s)
        # stopping condition = no more reactions in graph
        if len(_rxn_graph.reactions) == 0:
            return [rxn_graph.species[rxn_graph._specie_ids[s]] for s in list(mvc)]


if __name__ == "__main__":
    # reading graph from file
    rxn_graph = get_graph(settings.rxn_graph_path, to_networkx=False)
    print("Total number of species:", len(rxn_graph.species))
    print("Total number of reactions:", len(rxn_graph.reactions))
    G = rxn_graph.to_networkx_graph()
    # getting mvc via networkx
    networkx_mvc_smiles = networkx_mvc_algorithm(G)
    networkx_mvc_df = distance_distribution(G, [s.identifier for s in rxn_graph.reactant_species], networkx_mvc_smiles)
    networkx_mvc_df.to_csv(os.path.join(settings.parent_results_dir, "mvc_algs", "networkx_mvc.csv"))
    # getting mvc with greedy methods
    metrics = {"degree": total_degree, "percolation": percolation_degree}
    for name, metric in metrics.items():
        greedy_mvc_smiles = [s.identifier for s in greedy_mvc(rxn_graph, metric)]
        greedy_mvc_df = distance_distribution(G, [s.identifier for s in rxn_graph.reactant_species], greedy_mvc_smiles)
        greedy_mvc_df.to_csv(os.path.join(settings.parent_results_dir, "mvc_algs", "greedy_mvc_{}.csv".format(name)))
