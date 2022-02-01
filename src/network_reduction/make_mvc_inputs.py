# script to make initial geometry, given SMILES string of a specie and a target catalyst

import os
import openbabel as ob
import json
import numpy as np
import argparse
from .. import DirParser
from ..commons import get_graph, get_min_cover_species
import settings
from DirParser.src.utils.obUtils import obmol_to_molecule, ob_read_file_to_molecule, molecule_to_obmol
from DirParser.src.Base.Molecule import Molecule
from DirParser.src.Base.Bond import Bond


def build_molecule(catalyst: str, smiles_string: str):
    """Method to get the initial geometry guess for the molecule (specie or specie + catalyst).
    RETURNS (DirParser.Molecule) a DirParser molecule object with the guessed geometry"""
    # reading molecule to ob.OBMol
    obmol = ob.OBMol(); conv = ob.OBConversion()
    conv.SetInFormat("smi")
    conv.ReadString(obmol, smiles_string)
    # initial geometry guess (not affected by dummy atom)
    obmol = guess_geometry(obmol)
    # split build process for specie + catalyst and specie alone
    if "*" in smiles_string and len(smiles_string) == 1: # isolating case of "*" - only catalyst
        cat = get_catalyst(catalyst)
        return molecule_to_obmol(cat)
    elif "*" in smiles_string:
        # finding and removing dummy atom
        dummy_idx, active_atom_idx = get_dummy_atom_idx(obmol)
        obmol.DeleteAtom(obmol.GetAtom(dummy_idx + 1))
        if dummy_idx < active_atom_idx:
            active_atom_idx -= 1
        # joining specie with catalyst
        cat = get_catalyst(catalyst)
        mol = obmol_to_molecule(obmol)
        mol.properties["active_atom"] = active_atom_idx
        joined = join_molecules(cat, mol)
        return molecule_to_obmol(joined)
    else:
        return obmol

def get_dummy_atom_idx(obmol):
    """Method to mark the active atom in the specie and remove it. RETURNS dummy and active atom indicis (tuple)"""
    # checks if there is a dummy atom
    dummy = [atom.GetAtomicNum() for atom in ob.OBMolAtomIter(obmol)].index(0)
    active_atom = None
    for bond in ob.OBMolBondIter(obmol):
        if dummy == bond.GetBeginAtomIdx() - 1:
            active_atom = bond.GetEndAtom()
        elif dummy == bond.GetEndAtomIdx() - 1:
            active_atom = bond.GetBeginAtom()
    if not active_atom is None:
        return dummy, active_atom.GetIdx() - 1
    else:
        return dummy, 0

# define custom OpenbabelBuildError
class OpenbabelBuildError (Exception):
    pass

def guess_geometry(obmol):
    """Method to read the SMILES string to openbabel format and """
    # reading string to OBMol
    # generating geometry
    builder = ob.OBBuilder()
    if not builder.Build(obmol):
        conv = ob.OBConversion()
        raise OpenbabelBuildError("Failed building the smiles: ".format(conv.WriteString(obmol)))
    return obmol

def get_catalyst(catalyst: str):
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(parent_dir, "data", "catalysts_structures")
    with open(os.path.join(data_dir, "catalyst_info.json"), "r") as f:
        active_atom = json.load(f)[catalyst]
    cat = ob_read_file_to_molecule(os.path.join(data_dir, catalyst + ".xyz"))
    cat.properties["active_atom"] = active_atom
    return cat

def gen_translation_vector(catalyst: Molecule):
    """Method to make a normal vector to the 4N plane"""
    active = catalyst.atoms[catalyst.properties["active_atom"]]
    neighbors = catalyst.get_neighbors(catalyst.properties["active_atom"])
    v1 = active.coordinates - catalyst.atoms[neighbors[0]].coordinates
    v2 = active.coordinates - catalyst.atoms[neighbors[1]].coordinates
    res = np.array([v1[1] * v2[2] - v1[2] * v2[1],
                        v1[2] * v2[0] - v1[0] * v2[2],
                        v1[0] * v2[1] - v1[1] * v2[0]])
    return res / np.linalg.norm(res)


def join_molecules(catalyst: Molecule, specie: Molecule):
    """Method to join both catalyst and specie and make an initial guess for the geometry"""
    # centralizing both molecules
    catalyst.move_by_vector(np.negative(catalyst.atoms[catalyst.properties["active_atom"]].coordinates))
    translation_vec = np.negative(specie.atoms[specie.properties["active_atom"]].coordinates) + gen_translation_vector(catalyst)
    specie.move_by_vector(translation_vec)
    catalyst.join(specie)
    catalyst.add_bond(Bond(catalyst.properties["active_atom"], specie.properties["active_atom"] + len(catalyst.atoms) - len(specie.atoms), 1))
    return catalyst


def mm_geometry_optimization(obmol: ob.OBMol, force_field: str="UFF", n_steps: int=1000):
    '''Method to run a molecular mechanics geometry optimization. RETURNS openbabel molecule with optimized geometry.'''
    # optimization
    OBFF = ob.OBForceField.FindForceField(force_field)
    suc = OBFF.Setup(obmol)
    if not suc == True:
        raise RuntimeError("Could not set up force field for molecule")
    OBFF.ConjugateGradients(n_steps)
    OBFF.GetCoordinates(obmol)
    return obmol

def generate_geometry(catalyst: str, smiles_str: str):
    """Method to generate geometry (initial guess + MM optimization) for a given TorinaNet.Specie object. 
    ARGS:
        - catalyst (str): the catalyst name
        - smiles_str (str): smiles string of target specie to calculate
    RETURNS (DirParser.Molecule) a DirParser molecule object, to further run DFT calculation"""
    print(smiles_str.strip())
    # guessing specie geometry
    mol = build_molecule(catalyst, smiles_str)
    # optimizing geometry
    joined = mm_geometry_optimization(mol, "UFF", 10000)
    # returning DirParser.Molecule
    mol = obmol_to_molecule(joined)
    mol.properties["mm_energy"] = joined.GetEnergy(0)
    return mol

def print_msg(msg: str):
    l = max(int(round(len(msg) * 1.5)), 50)
    pre_space = int(round(l / 2 - len(msg) / 2))
    print("=" * l)
    print(" " * pre_space + msg)
    print("=" * l)


def make_orca_input(molecule: Molecule):
    pass

def test():
    rxn_graph_1 = get_graph(settings.rxn_graph_path, to_networkx=False)
    G1 = rxn_graph_1.to_netwokx_graph()
    reactants = [s.identifier for s in rxn_graph_1.reactant_species]
    mvc_species_1 = get_min_cover_species(G1)
    
    with open("./test.txt", "r") as f:
        old = []
        for line in f.readlines():
            old.append(line.strip())
    for s in old:
        if not s in mvc_species_1:
            print("DIFFERENT ! ({})".format(s))
    for s in mvc_species_1:
        if not s in old:
            print("DIFFERENT ! ({})".format(s))
    import sys; sys.exit()

    
    
if __name__ == "__main__":
    # making directories
    if not os.path.isdir(os.path.join(settings.parent_results_dir, "mvc", "structures")):
        os.makedirs(os.path.join(settings.parent_results_dir, "mvc", "structures"))
    # reading graph from file
    rxn_graph = get_graph(settings.rxn_graph_path, to_networkx=False)
    G = rxn_graph.to_netwokx_graph()
    reactants = [s.identifier for s in rxn_graph.reactant_species]
    # finding minimal vertex cover
    mvc_species = get_min_cover_species(G)
    # for every catalyst, make a guess structure for all MVC species
    for catalyst in ["cobalt-corrole"]:
        for specie in mvc_species:
            if not specie in reactants:
                mol = generate_geometry(catalyst, specie)
                name = "{}_{}.xyz".format(specie.replace("*", "M"), catalyst) if "*" in specie else "{}.xyz".format(specie)
                mol.save_to_file(os.path.join(settings.parent_results_dir, "mvc", "structures", name))
        for specie in reactants:
            mol = generate_geometry(catalyst, specie)
            mol.save_to_file(os.path.join(settings.parent_results_dir, "mvc", "structures", "{}.xyz".format(specie.replace("*", catalyst))))
