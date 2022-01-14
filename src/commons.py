import networkx as nx
from networkx.algorithms.approximation import vertex_cover
import numpy as np
import random
from __init__ import TorinaNet, DirParser
from TorinaNet.src.core.Specie import Specie as Specie
from TorinaNet.src.core.AcMatrix.BinaryAcMatrix import BinaryAcMatrix
from TorinaNet.src.core.RxnGraph import RxnGraph

def get_graph(path, to_networkx=True):
    read_graph = RxnGraph()
    read_graph.from_file(path)
    if to_networkx:
        return read_graph.to_netwokx_graph()
    else:
        read_graph.set_reactant_species([Specie.from_ac_matrix(BinaryAcMatrix.from_specie(Specie("[O][C][O]"))), 
                                            Specie.from_ac_matrix(BinaryAcMatrix.from_specie(Specie("[H]"))),
                                            Specie.from_ac_matrix(BinaryAcMatrix.from_specie(Specie("*")))])
        return read_graph


def project_to_type(G, target_type):
    new_G = nx.Graph()
    for node in G:
        if type(node) is target_type:
            neighbors = [nn for n in nx.all_neighbors(G, node) for nn in nx.all_neighbors(G, n)]
            for n in neighbors:
                new_G.add_edge(node, n)
    return new_G

def get_min_cover_species(G, add_weight=False):
    random.seed(0)
    np.random.seed(0)
    nG = project_to_type(G, str)
    if add_weight:
        nx.classes.function.set_node_attributes(nG, {n: 100 if "*" in n else 0 for n in nG}, "cost")
        ws = "cost"
    else:
        ws = None
    return vertex_cover.min_weighted_vertex_cover(nG, ws)
