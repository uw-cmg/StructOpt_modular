import random
import numpy as np
from structopt.tools import NeighborList

def rich2poor(atoms, max_natoms=0.05, surf_CN=11, factor=1.1):
    '''Moves an atom from a rich region to a poor region'''

    if type(max_natoms) is float:
        max_natoms = int(max_natoms * len(atoms))

    moves = random.randint(0, max_natoms)

    syms = np.asarray(atoms.get_chemical_symbols())
    NN_list = NeighborList(atoms, cutoff=None, factor=1.1)
    NNs = [list(syms[i]) for i in NN_list]
    unique_syms = list(set(atoms.get_chemical_symbols()))

    for _ in range(moves):

        # Pick an atom in a rich environment
        count_syms = np.zeros([len(unique_syms), len(NNs)])
        for i, sym in enumerate(unique_syms):
            count_syms[i] = [neighbors.count(sym)/len(neighbors) for neighbors in NNs]
        max_counts = np.max(count_syms.T, axis=1)
        prob_rich = max_counts / sum(max_counts)
        index_rich = np.random.choice(np.arange(len(syms)), p=prob_rich)
        element_rich = atoms[index_rich].symbol

        # Pick an atom of different species in poor environment of the same atom
        indices_poor = [i for i, sym in enumerate(syms) if sym != element_rich]
        min_counts = count_syms[unique_syms.index(element_rich)]
        min_counts = np.asarray([count for count, sym in zip(min_counts, syms)
                                 if sym != element_rich])
        min_counts = np.absolute(min_counts - 1)
        prob_poor = min_counts / sum(min_counts)
        index_poor = np.random.choice(indices_poor, p=prob_poor)
        element_poor = atoms[index_poor].symbol

        # Update symbols and nearest neighborlist
        syms[index_rich] = element_poor
        syms[index_poor] = element_rich

        for i, neighbors in enumerate(NN_list):
            if index_rich in neighbors:
                NNs[i][list(neighbors).index(index_rich)] = element_poor
            if index_poor in neighbors:
                NNs[i][list(neighbors).index(index_poor)] = element_rich

    atoms.set_chemical_symbols(syms)

    return
