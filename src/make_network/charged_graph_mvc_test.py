#!/home/shaharpit/miniconda3/envs/general/bin/python
#SBATCH -n 8
#SBATCH -o ./output_max_itr_2.txt
#SBATCH -e ./output_max_itr_2.txt
#SBATCH --mem=64G
import os
import json
# import sys; sys.path.append("/home/shaharpit/Personal/CO2NetProject/src")
from .. import TorinaNet as tn
tn.iterate.kernel = "vanilla"
from tn.core import Specie
from src.commons import distance_distribution
from src.network_reduction.mvc_algrothms import greedy_mvc, total_degree, percolation_degree

def analyze_graph(rxn_graph):
    """Method to run test on single reaction graph"""
    print("Number of species in graph:", len(rxn_graph.species))
    print("Number of reactions in graph:", len(rxn_graph.reactions))
    res = {}
    res["n_reactions"] = len(rxn_graph.reactions)
    res["n_species"] = len(rxn_graph.species)
    G = rxn_graph.to_networkx_graph()
    # making distance from source distribution
    print("make distance distribution...")
    distances_df = distance_distribution(G, [s.identifier for s in rxn_graph.reactant_species])
    res["distance_distribution"] = distances_df.to_dict()
    print("making mvc analysis...")
    # making MVCs with different metrics
    metrics = {"degree": total_degree, "percolation": percolation_degree}
    for name, metric in metrics.items():
        greedy_mvc_smiles = [s.identifier for s in greedy_mvc(rxn_graph, metric)]
        greedy_mvc_df = distance_distribution(G, [s.identifier for s in rxn_graph.reactant_species], greedy_mvc_smiles)
        res["{}_mvc_distribution".format(name)] = greedy_mvc_df.to_dict()
    print("done mvc analysis")
    print("ANALYSIS DONE")
    return res

max_bonds_dict = {
    0: 1,
    1: 1,
    6: 4,
    8: 2,
}

def max_bonds_per_atom(ac_matrix):
    """Filter that makes sure that atoms don't have too many bonds in the ac matrix"""
    for i in range(len(ac_matrix)):
        counter = 0
        atom = ac_matrix.get_atom(i)
        for j, bond in enumerate(ac_matrix.matrix[i]):
            if not j == i:
                # do not count bonds with catalyst for hydrogen and oxygen
                if (atom == 1 or atom == 8):
                    if not ac_matrix.get_atom(j) == 0:
                        counter += bond
                else:
                    counter += bond
        if max_bonds_dict[ac_matrix.get_atom(i)] < counter:
            return False
    return True

max_atoms_dict = {
    0: 1,
    4: 4,
    8: 4
}

def max_atoms(ac_matrix):
    counter = {0: 0, 4: 0, 8: 0}
    for i in range(len(ac_matrix)):
        atom = ac_matrix.get_atom(i)
        if not atom in counter:
            pass
        else:
            counter[atom] += 1
    return all([counter[k] <= max_atoms_dict[k] for k in max_atoms_dict.keys()])
    
def make_network(reactants, max_reduction, max_oxidation, max_abs_charge, enumerate_charges):
    rxn_graph = tn.core.RxnGraph()
    source_species = []
    for r in reactants:
        s = rxn_graph.add_specie(r)
        source_species.append(s)
    rxn_graph.set_source_species(source_species)
    ac_filters = [max_bonds_per_atom, max_atoms]
    conversion_filters = [tn.iterate.conversion_matrix_filters.MaxChangingBonds(3), 
                            tn.iterate.conversion_matrix_filters.OnlySingleBonds()]
    stop_cond = tn.iterate.stop_conditions.MaxIterNumber(2)
    iterator = tn.iterate.Iterator(rxn_graph)
    rxn_graph = iterator.enumerate_reactions(conversion_filters, ac_filters, stop_cond, verbose=1)
    if enumerate_charges:
        charge_iterator = tn.iterate.ChargeIterator(rxn_graph, type(rxn_graph))
        return charge_iterator.enumerate_charges(max_reduction, 
                                                    max_oxidation, 
                                                    [tn.iterate.charge_filters.MaxAbsCharge(max_abs_charge)])
    else:
        return rxn_graph


if __name__ == "__main__":
    reactions_dict = {
        "co2_reduction": {
                            "reactants": [Specie("O=C=O", charge=0), Specie("[H]", charge=1), Specie("*", charge=0)], 
                            "max_reduction": 1, 
                            "max_oxidation": 0, 
                            "max_abs_charge": 1, 
                            "enumerate_charges": True
                        },
        "water_splitting": {
                            "reactants": [Specie("O", charge=0), Specie("*", charge=0)], 
                            "max_reduction": 0, 
                            "max_oxidation": 1, 
                            "max_abs_charge": 1, 
                            "enumerate_charges": True
                        },
        "ethane_comb": {
                            "reactants": [Specie("O=O"), Specie("CC")], 
                            "max_reduction": None, 
                            "max_oxidation": None, 
                            "max_abs_charge": None, 
                            "enumerate_charges": False
                        },
        "propane_comb": {
                            "reactants": [Specie("O=O"), Specie("CCC")], 
                            "max_reduction": None, 
                            "max_oxidation": None, 
                            "max_abs_charge": None, 
                            "enumerate_charges": False
                        },
        "methane_comb": {
                            "reactants": [Specie("O=O"), Specie("C")], 
                            "max_reduction": None, 
                            "max_oxidation": None, 
                            "max_abs_charge": None, 
                            "enumerate_charges": False
                        },
        "ethene_comb": {
                            "reactants": [Specie("O=O"), Specie("C=C")], 
                            "max_reduction": None, 
                            "max_oxidation": None, 
                            "max_abs_charge": None, 
                            "enumerate_charges": False
                        }
    }

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("reaction", type=str, help="desired reaction to iterate ({})".format(", ".format(reactions_dict.keys())))
    reaction = parser.parse_args().reaction
    # enumerating reactions
    net = make_network(**reactions_dict[reaction.lower()])
    results_dir = "../results/mvc_test_13032022"
    with open(os.path.join(results_dir, "{}.json".format(reaction.lower()))) as f:
        json.dump(analyze_graph(net), f)
