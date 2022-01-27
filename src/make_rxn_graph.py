from datetime import date, datetime
from __init__  import TorinaNet
from TorinaNet.src.Iterate.daskIterator import Iterator
from TorinaNet.src.Iterate.conversion_matrix_filters import MaxChangingBonds, OnlySingleBonds
from TorinaNet.src.core.Specie import Specie

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
        for j, atom in enumerate(ac_matrix.matrix[i]):
            if not j == i:
                # do not count bonds with catalyst for hydrogen and oxygen
                if not ((atom == 1 or atom == 8) and ac_matrix.get_atom(i) == 0):
                    counter += 1
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
    


if __name__ == "__main__":
    reactants = [Specie("O=C=O"), Specie("[H]"), Specie("*")]
    max_itr = 2
    ac_filters = [max_bonds_per_atom, max_atoms]
    conversion_filters = [MaxChangingBonds(4), OnlySingleBonds()]
    print("Starting to iterate bonds with max {} macro-iterations".format(max_itr))
    iterator = Iterator()
    rxn_graph = iterator.gererate_rxn_network(reactants, max_itr, ac_filters, conversion_filters, verbose=1)
    rxn_graph.save_to_file("../data/rxn_net_{}.rxn".format(datetime.now().strftime("%d%m%Y")))
