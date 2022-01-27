from numpy import source
from commons import get_graph, distance_from_source

if __name__ == "__main__":
    paths = ["../data/rxn_net_15012022.rxn", "../data/rxn_net_16012022.rxn"]
    first_layer_species = []
    source_dist_tables = []
    for path in paths:
        print("running with", path)
        rxn_graph = get_graph(path, False)
        G = rxn_graph.to_networkx_graph()
        source_dists = distance_from_source(G, [s.identifier for s in rxn_graph.reactant_species])
        source_dist_tables.append(source_dists)
        first_layer_species.append(source_dists[source_dists["dist"] == 1].index)

    species1 = first_layer_species[0]
    species2 = first_layer_species[1]
    for s in species2:
        if s not in species1:
            print(s, source_dist_tables[0].loc[s, "dist"])