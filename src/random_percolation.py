# script to run random percolation test on the whole reaction graph
import os
from commons import percolation_test, get_graph
import settings

if __name__ == "__main__":
    rxn_graph = get_graph(settings.rxn_graph_path, to_networkx=False)
    print("Total number of species:", len(rxn_graph.species))
    print("Total number of reactions:", len(rxn_graph.reactions))
    res = percolation_test(rxn_graph, range(0, 100, 10), 5, None)
    res.to_csv(os.path.join(settings.parent_results_dir, "property_distributions", "population_percolation.csv"))