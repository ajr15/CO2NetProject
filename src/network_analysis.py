# script to visualize the network's species property distributions and the distance from the source
import os
from matplotlib import pyplot as plt
import networkx as nx
import pandas as pd
from copy import copy
from commons import get_graph, get_min_cover_species
import settings
from TorinaNet.src.analyze.visualize import specie_func_distribution
from TorinaNet.src.core.Specie import Specie



def calc_molar_mass(specie: Specie):
    obmol = specie.parse_identifier()
    return obmol.GetMolWt()

def calc_number_of_atoms(specie: Specie):
    obmol = specie.parse_identifier()
    return obmol.NumAtoms()

def shortest_path(G: nx.DiGraph, source, target):
    """Finds the shortest path between two nodes (species or reactions). The shortest path is given in terms of reaction nodes."""
    path = [x for x in nx.shortest_path(G, source, target) if type(x) is int]
    npaths = [path]
    rxns = []
    while True:
        # adding the reactions from the shortest path (greedy algorithm)
        nrxns = [i for i in npaths[npaths.index(min(npaths))] if type(i) is int]
        new_reactions = False
        for rxn in nrxns:
            if not rxn in rxns:
                rxns.append(rxn)
                new_reactions = True
        # stopping condition = no new reactions
        if not new_reactions:
            break
        # finding shortest paths for the new reactions (including longer creation paths for reactants)
        npaths = []
        for r in nrxns:
            for s in [i for i in G.neighbors(r) if not i == source or not i == target]:
                npaths.append([x for x in nx.shortest_path(G, source, s) if type(x) is int])
    return rxns


def distance_distribution(G, reactants, species=None):
    d = {}
    if species is None:
        species = [node for node in G if type(node) is str and not node in reactants]
    for specie in species:
        dist = min([len(shortest_path(G, r, specie)) for r in reactants])
        if not dist in d:
            d[dist] = 1
        else:
            d[dist] += 1
    return pd.DataFrame({"distantce": d.keys(), "count": d.values()})



if __name__ == "__main__":
    # reading graph from file
    rxn_graph = get_graph(settings.rxn_graph_path, to_networkx=False)
    G = rxn_graph.to_netwokx_graph()
    # analyzing shortest paths
    population_df = distance_distribution(G, [s.identifier for s in rxn_graph.reactant_species])
    population_df.to_csv(os.path.join(settings.parent_results_dir, "property_distributions", "distance_disribution.csv"))
    mvc_df = distance_distribution(G, [s.identifier for s in rxn_graph.reactant_species], get_min_cover_species(G))
    mvc_df.to_csv(os.path.join(settings.parent_results_dir, "property_distributions", "mvc_distance_disribution.csv"))
    # plot specie molar mass distribution
    specie_func_distribution(rxn_graph, calc_molar_mass, title="Molar Mass Distribution", x_label="Molar mass", njobs=4, color="black")
    plt.savefig(os.path.join(settings.parent_results_dir, "property_distributions", "molar_mass.png"))
    plt.close()
    # plot specie number of atoms distribution
    specie_func_distribution(rxn_graph, calc_number_of_atoms, title="Number of Atoms Distribution", x_label="Number of atoms", njobs=4, color="black")
    plt.savefig(os.path.join(settings.parent_results_dir, "property_distributions", "number_of_atoms.png"))
    plt.close()

