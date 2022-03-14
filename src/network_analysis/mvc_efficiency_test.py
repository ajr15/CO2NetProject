# script to test the added value of the MVC
# tested on generated networks of alkanes of various sizes
# 2 macro-iterations were used for all networks
# script used for network generation is in ../make_network/rxn_graphs_for_mvc_test.py
import os
import pandas as pd
from matplotlib import pyplot as plt
import sys; sys.path.append(".")
from src.commons import get_graph, distance_distribution
from src.network_reduction.mvc_algrothms import greedy_mvc, total_degree, percolation_degree

def test_mvc(rxn_graph, res_dir, prefix):
    """Method to run test on single reaction graph"""
    print("Running computation for {}".format(prefix))
    print("Number of species in graph:", len(rxn_graph.species))
    print("Number of reactions in graph:", len(rxn_graph.reactions))
    G = rxn_graph.to_networkx_graph()
    # making distance from source distribution
    print("make distance distribution...")
    distances_df = distance_distribution(G, [s.identifier for s in rxn_graph.reactant_species])
    distances_df.to_csv(os.path.join(res_dir, "{}_distance_dist.csv".format(prefix)))
    print("done distance distribution, saved at {}".format(os.path.join(res_dir, "{}_distance_dist.csv".format(prefix))))
    print("making mvc analysis...")
    # making MVCs with different metrics
    metrics = {"degree": total_degree, "percolation": percolation_degree}
    for name, metric in metrics.items():
        greedy_mvc_smiles = [s.identifier for s in greedy_mvc(rxn_graph, metric)]
        greedy_mvc_df = distance_distribution(G, [s.identifier for s in rxn_graph.reactant_species], greedy_mvc_smiles)
        greedy_mvc_df.to_csv(os.path.join(res_dir, "{}_greedy_mvc_{}.csv".format(prefix, name)))
    print("done mvc analysis")
    print("COMPUTATION DONE")

def analyze_results(res_dir):
    # reading data
    solution_sizes = {
        "naive": {},
        "degree_mvc": {},
        "percolation_mvc": {}
    }
    for res_file in os.listdir(res_dir):
        name = res_file.split("_")[0]
        df = pd.read_csv(os.path.join(res_dir, res_file))
        if "mvc_degree" in res_file:
            solution_sizes["degree_mvc"][name] = sum(df["count"].values)
        elif "mvc_percolation" in res_file:
            solution_sizes["percolation_mvc"][name] = sum(df["count"].values)
        elif "distance_dist" in res_file:
            solution_sizes["naive"][name] = df[df["distantce"] == 1]["count"]
    # plotting 
    species = sorted(solution_sizes["naive"].keys(), key=lambda x: len(x))
    for comp in solution_sizes.keys():
        plt.plot(range(1, len(species) + 1), [solution_sizes[comp][k] for k in species], label=comp)
    plt.legend()
    plt.title("MVC vs Naive Cover Sizes")
    plt.xlabel("# of carbons")
    plt.ylabel("# of species")
    plt.show()


    
if __name__ == "__main__":
    analyze_results("results/mvc_test_09022022/")
    sys.exit()
    rxn_graph_files = os.listdir("data/mvc_test/")
    for rxn_file in rxn_graph_files:
        reactant_smiles = rxn_file.split("_")[0]
        rxn_graph = get_graph(os.path.join("data/mvc_test/", rxn_file), False, [reactant_smiles])
        test_mvc(rxn_graph, "results/mvc_test_09022022/", reactant_smiles)