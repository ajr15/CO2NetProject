import networkx as nx
import numpy as np
import pandas as pd
from copy import copy

def distance_distribution(G, reactants, species=None):
    d = {}
    if species is None:
        species = [node for node in G if type(node) is str and not node in reactants]
    for specie in species:
        dists = []
        for r in reactants:
            path = shortest_path(G, r, specie)
            if not path is None:
                dists.append(len(path))
        dist = min(dists)
        d[specie] = dist
    return pd.DataFrame({"distantce": d.keys(), "count": d.values()})

def _get_path_reactants(G, target, path):
    """Method to get the reactants of the reactions in a path"""
    reactants = set()
    for rxn in path:
        for specie in G.neighbors(rxn):
            reactants.add(specie) 
    return reactants

def nx_shortest_path(G: nx.DiGraph, source, target, source_species):
    """Finds the shortest path between two nodes (species or reactions). The shortest path is given in terms of reaction nodes."""
    try:
        path = [x for x in nx.shortest_path(G, source, target) if type(x) is int]
    except nx.NetworkXNoPath:
        return None
    npaths = [path]
    rxns = []
    reactants = set()
    print("running for {} to {}".format(source, target))
    while True:
        # adding the reactions from the shortest path (greedy algorithm)
        nrxns = [i for i in npaths[npaths.index(min(npaths))] if type(i) is int]
        reactants = reactants.union(_get_path_reactants(nrxns))
        for rxn in nrxns:
            if not rxn in rxns:
                rxns.append(rxn)
        # stopping condition = no new reactions
        if all([a in source_species for a in list(reactants)]):
            break
        # finding shortest paths for the new reactions (including longer creation paths for reactants)
        npaths = []
        for r in nrxns:
            for s in list(reactants):
                try:
                    path = [x for x in nx.shortest_path(G, source, s) if type(x) is int]
                except nx.NetworkXNoPath:
                    continue
                npaths.append(path)
    return rxns

def distance_from_source(G, source_species):
    """Method to implement modified dijkstra's algorithm for finding distance between specie to the source of the network (reactants)"""
    # initializing 
    specie_df = {s: {"dist": np.inf, "rxn": None, "visited": False} for s in G if type(s) is str} # dictionary for each specie its distance from source and making reaction
    for s in source_species:
        specie_df[s] = {"dist": 0, "rxn": None, "visited": False}
    specie_df = pd.DataFrame.from_dict(specie_df, orient="index")
    # running main loop
    while not all(specie_df["visited"].values): # runs until all species are visited
        # find the next specie to visit (unvisited specie with shortest distance from source)
        unvisited = specie_df[specie_df["visited"] == False]
        specie = unvisited[unvisited["dist"] == unvisited["dist"].min()].first_valid_index()
        # chaning flag to visited
        specie_df.loc[specie, "visited"] = True
        # go over all reactions with visited reactants
        for rxn in G.successors(specie):
            pred_species = [s for s in G.predecessors(rxn)]
            if all([specie_df.loc[s, "visited"] for s in pred_species]):
                # make a distance estimate for the reaction's products
                # the estimate is the maximal distance of the reactant species plus one
                prod_species = [s for s in G.successors(rxn)]
                dist_estimate = max([specie_df.loc[s, "dist"] for s in pred_species]) + 1
                for s in prod_species:
                    # if the distance estimate is less than the known distance, update distance and pred reaction
                    if dist_estimate < specie_df.loc[s, "dist"]:
                        specie_df.loc[s, "dist"] = dist_estimate
                        specie_df.loc[s, "rxn"] = rxn
    return specie_df


def test_shortest_path():
    G = nx.DiGraph()
    G.add_edges_from([("source1", 1), 
                        (1, "mid1"),
                        ("mid1", 4),
                        ("source1", 4),
                        (4, "source2"),
                        (1, "mid2"),
                        ("mid1", 2),
                        ("source3", 2),
                        (2, "mid3"),
                        ("source2", 3), 
                        (3, "mid1"),
                        ("mid2", 5),
                        ("mid3", 5),
                        (5, "mid4"),
                        ("source3", 6),
                        ("mid3", 6),
                        (6, "mid2")])
    res = distance_from_source(G, ["source1", "source2", "source3"])
    print(res)
    import sys; sys.exit()

if __name__ == "__main__":
    test_shortest_path()
