"""
@uthor: Himaghna, 3rd December 2019
Description: Generate subgraphs from smiles

"""
from argparse import ArgumentParser
from os import getcwd, mkdir
import os.path
import pickle
from time import time

import networkx as nx
import networkx.algorithms.isomorphism as iso
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdmolops

def get_graph_descriptors(smiles, min_graph_length, max_graph_length):
    """
    Generate graph descriptors from smile strings

    Params ::
    smiles: np ndarray: collection of smiles of molecules
    min_graph_length: int: number of edges in the smallest sub-graph mined
    max_graph_length: int: number of edges in the largest sub-graph mined

    Returns ::
    dict with keys
        X: n x n numpy ndarray: Data matrix
        graph _descriptors: n x 1 List: list of sub-graph objects
            corresponding to columns of X
        failed_to_load_idx: list: indexes for which smiles could not be loaded
    """
    X, descriptors = list(), list()
    start_time = time()
    # some methods needed to generate graph descriptors
    def bond2graph(mol, bondids):

            """
            :param mol: RDKIT molecule object
            :param bonidx: tuple of bondids making the graph
            :return: corresponding networkx graph
            """
            g = nx.Graph()
            for bondid in bondids:
                Bond = mol.GetBondWithIdx(bondid)
                type_of_bond = Bond.GetBondType()
                begin_atom_idx = Bond.GetBeginAtomIdx()
                end_atom_idx = Bond.GetEndAtomIdx()
                g.add_node(begin_atom_idx,
                            atomic_number=(mol.GetAtomWithIdx(begin_atom_idx))
                                    .GetAtomicNum(),
                            atomic_symbol=(mol.GetAtomWithIdx(begin_atom_idx))
                                    .GetSymbol())
                g.add_node(end_atom_idx,
                            atomic_number=(mol.GetAtomWithIdx(end_atom_idx))
                                    .GetAtomicNum(),
                            atomic_symbol=(mol.GetAtomWithIdx(end_atom_idx))
                                    .GetSymbol())
                g.add_edge(begin_atom_idx, end_atom_idx, type=type_of_bond)

            return g

    def check_iso(graph1, graph2):
        """
        Check isomorphism using node attribute {atomic_number} and
        edge attribute {type}
        :param graph1: networkx graph Object
        :param graph2: networkx graph Object
        :return:  True is graph1 and graph 2 are isomorphic
        """
        return nx.is_isomorphic(graph1, graph2,
                                node_match=iso.categorical_node_match(
                                    'atomic_number', 0),
                                edge_match=iso.categorical_edge_match('type',
                                                                    'Single'))

    def update_X_descriptors(X, descriptors, graphs):
        """
        Updates coefficient matrix with counts of graph descriptors and creates
        new column for descriptor not seen before
        :param X: (2D list) Coefficient matrix with each row for  a molecule and
                each column for a coefficient of descriptor in molecule
        :param descriptor: (list) Unique graphs which form basis set of molecule
                        expansion
        :param graphs: (list) sub-graphs of a molecule
        :return: (updated) X, descriptors
        """
        row_entry = list()  # row which will become entry in X for this molecule
        row_entry.extend([0] * len(descriptors))
        for graph in (graphs):
            graph_classified = False # if graph has been assigned a descriptor
            for index, descriptor in enumerate(descriptors):
                if not len(descriptor.nodes()) == len(graph.nodes()) or \
                        not len(descriptor.edges()) == len(graph.edges()):
                    continue
                if check_iso(graph, descriptor):
                    row_entry[index] +=1
                    graph_classified = True
                    break  # stop going thru descriptors
            if not graph_classified:
                # graph is a descriptor not seen previously
                descriptors.append(graph)
                row_entry.append(1)
        X.append(row_entry)
        return X, descriptors

    load_fail_idx = []
    for count, smile in enumerate(smiles, 1):
#        print('Processing {} ({}/{})'.format(smile, count, len(smiles)))

        mol = Chem.MolFromSmiles(smile)
        if not mol:
            print(f'{smile} could not be loaded')
            load_fail_idx.append(count-1)
            continue
        # sanitize
        rdmolops.Kekulize(mol)
        graphs = list()
        subgraphs = rdmolops.FindAllSubgraphsOfLengthMToN(mol,
                                                          min=min_graph_length,
                                                          max=max_graph_length,
                                                          useHs=False)


        for subgraph_degree in subgraphs:
            for subgraph in subgraph_degree:
               graphs.append(bond2graph(mol=mol, bondids=subgraph))
        X, descriptors = update_X_descriptors(X=X, descriptors=descriptors,
                                              graphs=graphs)
    assert len(X[-1]) == max([len(row) for row in X])

    # coefficients of missing descriptors set to zero
    def complete_rows(n):
        """
        Given a 2D array, make it a square, i.e. make all rows of length n
        by padding the missing entries with zero
        :return: (2D List) Squared Matrix
        """
        X_out = []
        for row in X:
            row.extend([0] * (n - len(row)))
            # row_entry.extend([0] * len(descriptors))
            X_out.append(row)
        return X_out

    X = complete_rows(n=len(X[-1]))
#    print(f'Time for max graph {max_graph_length} is {time() - start_time}')
    return {
        'X': np.array(X),
        'graph_descriptors': descriptors,
        'failed_to_load_idx': load_fail_idx
    }




if __name__ == "__main__":
    main()

