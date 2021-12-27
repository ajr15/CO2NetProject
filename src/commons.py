import os
import openbabel as ob
import json
import sys
import numpy as np
sys.path.append("C:\\Users\\Shachar\\OneDrive - Technion\\Technion\\PhD\\Complex Networks")
from TorinaNet.src.core.Specie import Specie as trSpecie
sys.path.append("C:\\Users\\Shachar\\OneDrive - Technion\\Technion\\PhD\\Other Projects")
from DirParser.src.utils.obUtils import obmol_to_molecule, ob_read_file_to_molecule, molecule_to_obmol
from DirParser.src.Base.Molecule import Molecule
from DirParser.src.Base.Bond import Bond


def guess_specie_geometry(specie: trSpecie):
    """Method to get the initial guess for the specie (organic part only).
    RETURNS (DirParser.Molecule) a DirParser molecule object with the guessed geometry"""
    # reading string to OBMol
    obmol = ob.OBMol(); conv = ob.OBConversion()
    conv.SetInFormat("smi")
    conv.ReadString(obmol, specie.identifier)
    # generating geometry
    builder = ob.OBBuilder()
    builder.Build(obmol)
    # find binding atom
    dummy = [atom.GetAtomicNum() for atom in ob.OBMolAtomIter(obmol)].index(0)
    for bond in ob.OBMolBondIter(obmol):
        if dummy == bond.GetBeginAtomIdx() - 1:
            active_atom = bond.GetEndAtom()
        elif dummy == bond.GetEndAtomIdx() - 1:
            active_atom = bond.GetBeginAtom()
    obmol.DeleteAtom(obmol.GetAtom(dummy + 1))
    # convert to DirParser.Molecule    
    mol = obmol_to_molecule(obmol)
    mol.properties["active_atom"] = active_atom.GetIdx() - 1
    print(mol.atoms[mol.properties["active_atom"]], mol.properties["active_atom"])
    return mol


def get_catalyst(catalyst: str):
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(parent_dir, "data", "catalysts_structures")
    with open(os.path.join(data_dir, "catalyst_info.json"), "r") as f:
        active_atom = json.load(f)[catalyst]
    cat = ob_read_file_to_molecule(os.path.join(data_dir, catalyst + ".xyz"))
    cat.properties["active_atom"] = active_atom
    return cat


def join_ob_mols(mol1, mol2):
    for atom in ob.OBMolAtomIter(mol2):
        mol1.AddAtom(atom)
    for bond in ob.OBMolBondIter(mol2):
        trail = mol1.AddBond(bond)
        if not trail:
            trail2 = mol1.AddBond(bond.GetBeginAtomIdx() + mol1.NumAtoms() - mol2.NumAtoms(), bond.GetEndAtomIdx() + mol1.NumAtoms() - mol2.NumAtoms(), bond.GetBO())
    return mol1

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
    objoined = molecule_to_obmol(catalyst)
    conv = ob.OBConversion()
    conv.WriteFile(objoined, "./guess.mol")
    return objoined


def mm_geometry_optimization(obmol: ob.OBMol, ncatalyst_atoms: int, force_field: str="UFF", n_steps: int=1000):
    '''Method to run a molecular mechanics geometry optimization. RETURNS openbabel molecule with optimized geometry.'''
    # optimization
    OBFF = ob.OBForceField.FindForceField(force_field)
    constraints = ob.OBFFConstraints()
    # for i in range(ncatalyst_atoms):
    #     constraints.AddAtomConstraint(i)
    suc = OBFF.Setup(obmol, constraints)
    if not suc == True:
        raise RuntimeError("Could not set up force field for molecule")
    OBFF.ConjugateGradients(n_steps)
    OBFF.GetCoordinates(obmol)
    return obmol

def generate_geometry(catalyst: str, specie: trSpecie):
    """Method to generate geometry (initial guess + MM optimization) for a given TorinaNet.Specie object. 
    ARGS:
        - catalyst (str): the catalyst name
        - specie (TorinaNet.Specie): target specie to calculate
    RETURNS (DirParser.Molecule) a DirParser molecule object, to further run DFT calculation"""
    # guessing specie geometry
    mol = guess_specie_geometry(specie)
    mol.save_to_file("./mol.xyz")
    # reading catalyst
    cat = get_catalyst(catalyst)
    ncat_atoms = len(cat.atoms)
    # making initial guess
    print("joining")
    joined = join_molecules(cat, mol)
    # optimizing geometry
    print("optimizing")
    joined = mm_geometry_optimization(joined, ncat_atoms, "UFF", 10000)
    # returning DirParser.Molecule
    return obmol_to_molecule(joined)

if __name__ == "__main__":
    specie = trSpecie("*c1c([H])c([H])c([H])c([H])c1([H])")
    print("specie =", specie.identifier)
    catalyst = "cobalt-corrole"
    print("generating geometry")
    guess = generate_geometry(catalyst, specie)
    print("saving to file")
    guess.save_to_file("./test.xyz")
    print("done !")