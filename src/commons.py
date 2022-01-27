import networkx as nx
from networkx.algorithms.approximation import vertex_cover
import numpy as np
import pandas as pd
import random
from __init__ import TorinaNet, DirParser
from TorinaNet.src.core.Specie import Specie as Specie
from TorinaNet.src.core.AcMatrix.BinaryAcMatrix import BinaryAcMatrix
from TorinaNet.src.core.RxnGraph import RxnGraph

def get_graph(path, to_networkx=True):
    read_graph = RxnGraph()
    read_graph.from_file(path)
    if to_networkx:
        return read_graph.to_netwokx_graph()
    else:
        read_graph.set_reactant_species([Specie.from_ac_matrix(BinaryAcMatrix.from_specie(Specie("[O][C][O]"))), 
                                            Specie.from_ac_matrix(BinaryAcMatrix.from_specie(Specie("[H]"))),
                                            Specie.from_ac_matrix(BinaryAcMatrix.from_specie(Specie("*")))])
        return read_graph


def distance_from_source(G, source_species):
    # TODO: integrate this method into TorinaNet.analysis package !
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
    # dropping the "visited column"
    specie_df = specie_df.drop(columns=["visited"])
    return specie_df


def get_path_to_source(G, distance_table, target):
    rxn_path = [distance_table.loc[target, "rxn"]]
    while True:
        parent_rxn = rxn_path[-1]
        # finding the predecessor reactant (one with max distance)
        max_dist = 0
        parent_specie = None
        for reactant in G.predecessors(parent_rxn):
            dist = distance_table.loc[reactant, "dist"]
            if dist > max_dist:
                parent_specie = reactant
                max_dist = dist
        # stoping condition = the predecessor is a source specie
        if max_dist == 0:
            return rxn_path
        # adding the predecessor's reaction to the rxn_path
        new_rxn = distance_table.loc[parent_specie, "rxn"]
        rxn_path.append(new_rxn)


def calc_betweeness_centrality(G, source_species):
    # init
    distance_table = distance_from_source(G, source_species)
    counter_dict = {s: 0 for s in G if type(s) is str and not s in source_species}
    # iterating over species
    for specie in G:
        if type(specie) is str and not specie in source_species:
            # find shortest path to source for each
            rxn_path = get_path_to_source(G, distance_table, specie)
            species = [s for rxn in rxn_path for s in nx.all_neighbors(G, rxn)]
            # add species on shortest path to counter
            for s in species:
                if not s in source_species and not s == specie:
                    counter_dict[s] += 1
    return counter_dict

def distance_distribution(G, reactants, species=None):
    d = {}
    distances = distance_from_source(G, reactants)
    if species is None:
        species = [node for node in G if type(node) is str]
    for specie in species:
        dist = distances.loc[specie, "dist"]
        if not dist in d:
            d[dist] = 1
        else:
            d[dist] += 1
    return pd.DataFrame({"distantce": d.keys(), "count": d.values()})

def random_percolation(rxn_graph, nspecies, nreps, species_set=None):
    """Method to make random percolation test for rxn_graph.
    removing nspecies with nreps"""
    avg_species = 0
    avg_rxns = 0
    avg_removes = 0
    print("PERCOLATION TEST START !")
    if species_set is None:
        species_set = rxn_graph.species
    for rep in range(nreps):
        print("RUNNING REP {} OUT OF {}".format(rep + 1, nreps))
        nremoves = 0
        for _ in range(nspecies):
            specie = species_set[np.random.randint(0, len(species_set) - 1)]
            if rxn_graph.has_specie(specie):
                try:
                    rxn_graph = rxn_graph.remove_specie(specie)
                    nremoves += 1
                except KeyError:
                    print("ERROR IN SPECIE", specie.identifier)
                    continue
        avg_species += len(rxn_graph.species)
        avg_rxns += len(rxn_graph.reactions)
        avg_removes += nremoves
    print("PERCOLATION TEST DONE !")
    return avg_species / nreps, avg_rxns / nreps, avg_removes / nreps


def percolation_test(rxn_graph, n_species_range, nreps, species_set):
    res = pd.DataFrame()
    for i, nspecies in enumerate(n_species_range):
        print("**** Running percolation test with {} removals ({} out of {})".format(nspecies, i + 1, len(n_species_range)))
        try:    
            avg_species, avg_rxns, avg_removes = random_percolation(rxn_graph, nspecies, nreps, species_set)
            res = res.append({"n_species": avg_species, "n_reactions": avg_rxns, "n_removes": avg_removes}, ignore_index=True)
        except ValueError:
            print("ERROR IN PERCOLATION TEST")
            continue
    res = res.set_index("n_removes")
    return res