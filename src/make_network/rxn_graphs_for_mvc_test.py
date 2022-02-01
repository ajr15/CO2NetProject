#!/home/shaharpit/miniconda3/envs/general/bin/python
#SBATCH -n 8
#SBATCH -o ./output_max_itr_2.txt
#SBATCH -e ./output_max_itr_2.txt
#SBATCH --mem=64G

from datetime import datetime
# import sys; sys.path.append("/home/shaharpit/Personal/CO2NetProject/src")
from .. import TorinaNet
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
    
def make_network(reactant_smiles: str):
    reactants = [Specie(reactant_smiles)]
    max_itr = 2
    ac_filters = [max_bonds_per_atom, max_atoms]
    conversion_filters = [MaxChangingBonds(3), OnlySingleBonds()]
    iterator = Iterator()
    rxn_graph = iterator.gererate_rxn_network(reactants, max_itr, ac_filters, conversion_filters, verbose=1)
    rxn_graph.save_to_file("data/mvc_test/{}_{}.rxn".format(reactant_smiles, datetime.now().strftime("%d%m%Y")))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("reactant", type=str, help="smiles string of the reactant")
    args = parser.parse_args()
    make_network(args.reactant)
