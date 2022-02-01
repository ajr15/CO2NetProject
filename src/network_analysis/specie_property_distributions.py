# script to calculate distributions of specie's properties in the graph
import os
from matplotlib import pyplot as plt
from __init__ import TorinaNet
from TorinaNet.src.analyze.visualize import specie_func_distribution
from TorinaNet.src.core.Specie import Specie
import settings
from commons import get_graph, distance_distribution


def calc_molar_mass(specie: Specie):
    obmol = specie.parse_identifier()
    wt = obmol.GetMolWt()
    if "*" in specie.identifier:
        return wt + 300
    else:
        return wt

def calc_number_of_atoms(specie: Specie):
    obmol = specie.parse_identifier()
    wt = obmol.NumAtoms()
    if "*" in specie.identifier:
        return wt + 100
    else:
        return wt

if __name__ == "__main__":
    print("Reading graph from", settings.rxn_graph_path)
    rxn_graph = get_graph(settings.rxn_graph_path, to_networkx=False)
    print("number of species in graph", len(rxn_graph.species))
    print("number of reactions in graph", len(rxn_graph.reactions))
    G = rxn_graph.to_networkx_graph()
    # analyzing shortest paths
    population_df = distance_distribution(G, [s.identifier for s in rxn_graph.reactant_species])
    population_df.to_csv(os.path.join(settings.parent_results_dir, "property_distributions", "distance_disribution.csv"))

    print("Plotting molar mass distribution")
    # plot specie molar mass distribution
    specie_func_distribution(rxn_graph, calc_molar_mass, title="Molar Mass Distribution", x_label="Molar mass", njobs=4, color="black", bins=100)
    plt.savefig(os.path.join(settings.parent_results_dir, "property_distributions", "molar_mass.png"))
    plt.close()
    # plot specie number of atoms distribution
    print("Plotting number of atoms distribution")
    specie_func_distribution(rxn_graph, calc_number_of_atoms, title="Number of Atoms Distribution", x_label="Number of atoms", njobs=4, color="black", bins=100)
    plt.savefig(os.path.join(settings.parent_results_dir, "property_distributions", "number_of_atoms.png"))
    plt.close()
    print("ALL DONE !")
