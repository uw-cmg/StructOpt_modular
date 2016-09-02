import random
import numpy as np
from structopt.tools import NeighborList

def poor2rich(atoms, max_natoms=0.05, surf_CN=11, factor=1.1):
    '''Moves an atom from a rich region to a poor region'''

    if type(max_natoms) is float:
        max_natoms = int(max_natoms * len(atoms))

    moves = random.randint(0, max_natoms)

    syms = np.asarray(atoms.get_chemical_symbols())
    NN_list = NeighborList(atoms, cutoff=None, factor=1.1)
    NNs = [list(syms[i]) for i in NN_list]
    unique_syms = list(set(atoms.get_chemical_symbols()))

    for _ in range(moves):

        # Pick an atom in a poor environment
        count_syms = np.zeros([len(unique_syms), len(NNs)])
        for i, sym in enumerate(unique_syms):
            count_syms[i] = [neighbors.count(sym)/len(neighbors) for neighbors in NNs]
        min_counts = np.min(count_syms.T, axis=1)
        prob_poor = np.absolute(min_counts - 1)
        prob_poor = prob_poor / sum(prob_poor)
        index_poor = np.random.choice(np.arange(len(syms)), p=prob_poor)
        element_poor = atoms[index_poor].symbol

        # Pick an atom of different species in rich environment of the same atom
        indices_rich = [i for i, sym in enumerate(syms) if sym != element_poor]
        max_counts = count_syms[unique_syms.index(element_poor)]
        max_counts = np.asarray([count for count, sym in zip(max_counts, syms)
                                 if sym != element_poor])
        prob_rich = max_counts / sum(max_counts)
        index_rich = np.random.choice(indices_rich, p=prob_rich)
        element_rich = atoms[index_rich].symbol

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
