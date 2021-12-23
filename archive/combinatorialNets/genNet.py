import sys; sys.path.append("/home/shachar/Documents/")
from TorinaNet.src.Iterate.Iterator import Iterator
from TorinaNet.src.core.Specie import Specie
from TorinaNet.src.Iterate.ac_matrix_filters import max_bonds_per_atom
from TorinaNet.src.Iterate.conversion_matrix_filters import MaxChangingBonds, OnlySingleBonds
from TorinaNet.src.core.RxnGraph import RxnGraph

if __name__ == '__main__':
    reactants = [Specie("O=C=O"), Specie("[H][H]")]
    # reactants = [Specie("O")]
    max_itr = 2
    ac_filters = [max_bonds_per_atom]
    conversion_filters = [MaxChangingBonds(3), OnlySingleBonds()]
    iterator = Iterator()

    print("Iterating over species...")
    rxn_graph = iterator.gererate_rxn_network(reactants, max_itr, ac_filters, conversion_filters)
    print("Number of species in combinatorial graph:", len(rxn_graph.species))
    print("Number of reactions in combinatorial graph:", len(rxn_graph.reactions))
    print("saving graph to file")
    rxn_graph.save_to_file("./rxn_graph_max3.dat")
    print("reading graph from file")
    read_graph = RxnGraph()
    read_graph.from_file("./rxn_graph.dat")
    print("Number of species in read graph:", len(read_graph.species))
    print("Number of reactions in read graph:", len(read_graph.reactions))
    print("ALL DONE !")
    