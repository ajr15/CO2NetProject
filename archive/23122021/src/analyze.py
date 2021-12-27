import sys
from networkx.algorithms import shortest_paths
import networkx as nx
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import os
from TorinaNet.src.core.RxnGraph import RxnGraph
from TorinaNet.src.core.Specie import Specie
from TorinaNet.src.core.AcMatrix.BinaryAcMatrix import BinaryAcMatrix



def get_graph(path, to_networkx=True):
    read_graph = RxnGraph()
    read_graph.from_file(path)
    print("Total number of reactions =", len(read_graph.reactions))
    print("Total number of species =", len(read_graph.species))
    if to_networkx:
        return read_graph.to_netwokx_graph()
    else:
        read_graph.set_reactant_species([Specie.from_ac_matrix(BinaryAcMatrix.from_specie(Specie("O=C=O"))), 
                                            Specie.from_ac_matrix(BinaryAcMatrix.from_specie(Specie("[H]"))),
                                            Specie.from_ac_matrix(BinaryAcMatrix.from_specie(Specie("*")))])
        return read_graph

def _make_degree_dist_plot(n_species, degree_list, title, add_poiss, rxn_types, image_dir, prefix):
    plt.figure()
    plt.title(title)
    n, bins, patches = plt.hist(degree_list, label='count', bins=20)
    plt.ylabel("Count")
    plt.xlabel("Degree")
    if add_poiss:
        # plotting by avg degree
        l = np.mean(degree_list)
        f = lambda x: (n_species) * l**x * np.exp(-l) / np.math.factorial(x)
        x_vals = [round(bins[0]) + i for i in range(int(round(bins[-1] - bins[0]) + 1))]
        f_vals = [f(x) for x in x_vals]
        plt.plot(x_vals, f_vals, c='k', label='erdos-renyi')
        # plotting by rxn types
        if not rxn_types is None:
            nrxns = np.sum([i for i in rxn_types.values()])
            p = nrxns / (n_species * (n_species - 1)) # probability of reaction between 2 spcies
            l = n_species / nrxns * (rxn_types['11'] + 0.3333 * (rxn_types['12'] + rxn_types['21']) + 0.5 * rxn_types['22']) * p
            f = lambda x: (n_species) * l**x * np.exp(-l) / np.math.factorial(x)
            x_vals = [round(bins[0]) + i for i in range(int(round(bins[-1] - bins[0]) + 1))]
            f_vals = [f(x) for x in x_vals]
            plt.plot(x_vals, f_vals, c='r', label='rxn model')
        plt.legend()
    if not image_dir is None:
        plt.savefig(os.path.join(image_dir, prefix + '.png'))

def plot_degree_distribution(G, t='all', add_poiss=True, image_dir=None, prefix=''):
    if not os.path.isdir(image_dir):
        os.mkdir(image_dir)
    if t not in ['all', 'r_degrees', 'p_degrees', 'tot_degrees']:
        raise ValueError("Unrecognized option for \'t\' keyword, allowed plots types are: all, r_degrees, p_degrees and tot_degrees")
    r_degrees = []
    p_degrees = []
    tot_degrees = []
    counter = 0
    for node in G:
        if type(node) is str:
            counter += 1
            in_d = len([x for x in G.neighbors(node)])
            out_d = len(list(nx.all_neighbors(G, node))) - in_d
            r_degrees.append(in_d)
            p_degrees.append(out_d)
            tot_degrees.append(in_d + out_d)
    if add_poiss:
        rxn_types = get_rxn_types(G)
    else:
        rxn_types = None
    if t == 'all' or t == 'r_degrees':
        _make_degree_dist_plot(counter, r_degrees, "Reactant degrees", add_poiss, rxn_types, image_dir, prefix + "_reactants_dist")
    if t == 'all' or t == 'p_degrees':
        _make_degree_dist_plot(counter, p_degrees, "Product degrees", add_poiss, rxn_types, image_dir, prefix + "_products_dist")
    if t == 'all' or t == 'tot_degrees':
        _make_degree_dist_plot(counter, tot_degrees, "Total degrees", add_poiss, rxn_types, image_dir, prefix + "_degree_dist")

def get_rxn_types(G):
    rxn_types = {'11': 0,
                '12': 0,
                '21': 0,
                '22': 0,
                'other': 0}
    for node in G:
        if type(node) is int:
            in_d = len([x for x in G.neighbors(node)])
            out_d = len(list(nx.all_neighbors(G, node))) - in_d
            if in_d == 1 and out_d == 1:
                rxn_types['11'] += 1
            elif in_d == 1 and out_d == 2:
                rxn_types['12'] += 1
            elif in_d == 2 and out_d == 1:
                rxn_types['21'] += 1
            elif in_d == 2 and out_d == 2:
                rxn_types['22'] += 1
            else:
                rxn_types['other'] += 1
    return rxn_types

def plot_rxn_types(G, image_dir=None, prefix='rxn_types'):
    rxn_types = get_rxn_types(G)
    plt.figure()
    plt.title("Reaction types")
    plt.ylabel("Count")
    plt.annotate("Total number of reactions: {}".format(np.sum([i for i in rxn_types.values()])), xy=(0.55, 0.95), xycoords='axes fraction')
    plt.bar([i + 1 for i in range(5)], [i for i in rxn_types.values()], tick_label=[i for i in rxn_types.keys()])
    if not image_dir is None:
        plt.savefig(os.path.join(image_dir, prefix + '.png'))            

def make_text_pos(pos, xshift, yshift):
    for k, v in pos.items():
        pos[k] = v - np.array([xshift, yshift])
    return pos

def project_to_type(G, target_type):
    new_G = nx.Graph()
    for node in G:
        if type(node) is target_type:
            neighbors = [nn for n in nx.all_neighbors(G, node) for nn in nx.all_neighbors(G, n)]
            for n in neighbors:
                new_G.add_edge(node, n)
    return new_G


def visualize_rxn_graph(G, image_dir=None, prefix="Graph"):
    print("before:", len(G))
    oldG = G.copy()
    G.remove_nodes_from([node for node in G if type(node) is str and len(list(nx.all_neighbors(G, node))) < 10])
    G = project_to_type(G, str)
    print("after:", len(G))
    colormap = []
    sizemap = []
    labels = {}
    for node in G:
        if type(node) is str:
            if len(list(nx.all_neighbors(oldG, node))) < 200:
                colormap.append("gray")
                sizemap.append(50)
            else:
                colormap.append("red")
                sizemap.append(50)
                labels[node] = len(list(nx.all_neighbors(G, node)))
        else:
            colormap.append("green")
            sizemap.append(100)
    plt.figure(figsize=(20, 20))
    pos = nx.spring_layout(G, k=1)
    nx.draw_networkx_nodes(G, pos=pos, node_size=sizemap, node_color=colormap)
    nx.draw_networkx_edges(G, pos=pos, alpha=0.04)
    # nx.draw_networkx_labels(G, labels=labels, pos=make_text_pos(pos, 0, 0.06), font_color='k')
    for k, v in labels.items():
        print(k[:-1], v)
    if not image_dir is None:
        plt.savefig(os.path.join(image_dir, prefix + '.png'))

def _shortest_path_walker(G, path):
    ajr = []
    for n in path:
        if type(n) is int:
            reactants = [r for r in G.neighbors(n)]
            for r in reactants:
                ajr.append(r)
    return ajr

def shortest_path(G, source, target):
    path = nx.shortest_path(G, source, target)
    species = set([n for n in path if type(n) is str and not n == source and not n == target])
    npaths = [path]
    while True:
        nspecies = []
        for path in npaths:
            for specie in _shortest_path_walker(G, path):
                if not specie in species and not specie == source and not specie == target:
                    nspecies.append(specie)
                    species.add(specie)
        if len(nspecies) == 0:
            break
        npaths = []
        for s in nspecies:
            npaths.append(nx.shortest_path(G, source, s))
    return species

def shortest_path_rank(G, reactants):
    d = {}
    counter = 0
    for node in G:
        if type(node) is str and not node in [s for s in reactants]:
            counter += 1
            spaths = [list(shortest_path(G, s, node[:-2])) for s in reactants]
            spath = spaths[[len(p) for p in spaths].index(max([len(p) for p in spaths]))]
            for s in list(spath):
                if not s in d:
                    d[s] = 1
                else:
                    d[s] += 1
    print("Scanned {} species".format(counter))
    return pd.DataFrame({"specie": d.keys(), "rank": d.values()})

def get_min_cover_species(G, add_weight=False):
    nG = project_to_type(G, str)
    if add_weight:
        nx.classes.function.set_node_attributes(nG, {n: Specie(n[:-2]).parse_identifier().NumAtoms() for n in nG}, "num_atoms")
        print(list(nG.nodes(data="num_atoms"))[100])
        ws = "num_atoms"
    else:
        ws = None
    return nx.algorithms.approximation.vertex_cover.min_weighted_vertex_cover(nG, ws)

def analyze_covers(add_w):
    path = "../data/rxn_net_30112021.rxn"
    G = get_graph(path)
    cover = get_min_cover_species(G, add_weight=add_w)
    print("Number of speices in cover:", len(list(cover)))
    print("\nDistance from source of cover species:\n")
    count_d = {}
    for n in list(cover):
        p1 = shortest_path(G, "[H]\t\n", n)
        p2 = shortest_path(G, "[O][C][O]\t\n", n)
        p3 = shortest_path(G, "*\t\n", n)
        l = max(len(p1), len(p2), len(p3))
        if not l in count_d:
            count_d[l] = 1
        else:
            count_d[l] += 1
        print(l, n[:-2])
    print("\nDistribution of distances\n")
    print("\n".join(["{} = {} ({} %)".format(k, v, round(v / sum(count_d.values()) * 100, 2)) for k, v in count_d.items()]))

def plot_func_value_distribution(rxn_graph, func, plot_title, x_label, save_path, **plot_kwargs):
    """Method to make a histogram of function values. function takes a Specie from rxn_graph"""
    plt.figure()
    vals = []
    for specie in rxn_graph.species:
        vals.append(func(specie))
    plt.hist(vals, **plot_kwargs)
    plt.title(plot_title)
    plt.xlabel(x_label)
    plt.ylabel("Count")
    plt.savefig(save_path)

def plot_num_atoms_distribution(rxn_graph, image_dir):
    
    def f(specie):
        obmol = specie.parse_identifier()
        return obmol.NumAtoms()

    plot_func_value_distribution(rxn_graph, f, "Number of Atoms", "Number of atoms", os.path.join(image_dir, "num_atoms.png"))

def plot_molar_mass_distribution(rxn_graph, image_dir):
    
    def f(specie):
        obmol = specie.parse_identifier()
        return obmol.GetMolWt()

    plot_func_value_distribution(rxn_graph, f, "Molar Mass", "Molar mass", os.path.join(image_dir, "molar_mass.png"))

def plot_shortest_path_ranks(G, image_dir):
    plt.figure()
    ranks = shortest_path_rank(G, ["[H]", "[O][C][O]", "*"])["rank"]
    plt.hist(ranks)
    plt.title("Shortest Paths Count")
    plt.savefig(os.path.join(image_dir, "shortest_path_rank.png"))


def random_percolation(reaction_graph: RxnGraph, n_removals):
    """Method to return the number of species and reactions (as tuple) of reaction graph after n random removals"""
    species = np.random.choice(reaction_graph.species, n_removals)
    res = reaction_graph.remove_specie(species[0])
    for specie in species[1:]:
        if res.has_specie(specie):
            res = res.remove_specie(specie)
    return len(res.species), len(res.reactions)


def percolation_test(reaction_graph, max_n, n_reps, csv_path, step_size=5):
    res = pd.DataFrame()
    for n in range(1, max_n + 1, step_size):
        print("running for n = {}".format(n))
        avg_species = 0
        avg_rxns = 0
        for _ in range(n_reps):
            nspecies, nrxns = random_percolation(reaction_graph, n)
            avg_species += nspecies
            avg_rxns += nrxns
        res = res.append({"n_removes": n, "n_species": avg_species / n_reps, "n_reactions": avg_rxns / n_reps}, ignore_index=True)
        res.to_csv(csv_path)
    return res


def make_plots(date_str: str):
    path = "../data/rxn_net_30112021.rxn"
    image_dir = "../results/{}".format(date_str)
    G = get_graph(path)
    plot_degree_distribution(G, image_dir=image_dir, add_poiss=False)
    analyze_covers(False)
    plot_shortest_path_ranks(G, image_dir)
    rxn_graph = get_graph(path, False)
    plot_num_atoms_distribution(rxn_graph, image_dir)
    plot_molar_mass_distribution(rxn_graph, image_dir)
    plt.show()

def main():
    date_str = "15122021"
    path = "../data/rxn_net_30112021.rxn"
    G = get_graph(path, to_networkx=True)
    # res = percolation_test(G, 50, 5, "../results/15122021/percolation.csv", step_size=5)
    dist_d = {}
    for node in G:
        if type(node) is str:
            dist = max(len(shortest_path(G, "[H]\t\n", node)), 
                        len(shortest_path(G, "[O][C][O]\t\n", node)),
                        len(shortest_path(G, "*\t\n", node)))
            if not dist in dist_d:
                dist_d[dist] = 1
            else:
                dist_d[dist] += 1
    print("Distance distribution")
    for k, v in dist_d.items():
        print("{}: {}".format(k, v))
    analyze_covers(True)



if __name__ == '__main__':
    main()