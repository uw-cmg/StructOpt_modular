import random
import numpy as np
from copy import deepcopy

from structopt.tools import CoordinationNumbers

def swap_core_shell(individual, max_natoms=0.2, surf_CN=11):
    """Swaps atoms on the surface with an atom in the core. Only does it
    for different element types"""

    if not len(individual):
        return None

    CNs = CoordinationNumbers(individual)

    surf_indices = [i for i, CN in enumerate(CNs) if CN < surf_CN and CN > 3]
    bulk_indices = [i for i in range(len(individual)) if i not in surf_indices]
    print('n surface atoms:', len(surf_indices))
          
    # Construct surface and bulk dictionaries of elements and their indices
    syms = individual.get_chemical_symbols()
    surf_elements = list(set([syms[i] for i in surf_indices]))
    bulk_elements = list(set([syms[i] for i in bulk_indices]))

    surf_dict = {element: [] for element in surf_elements}
    bulk_dict = {element: [] for element in bulk_elements}

    for i in surf_indices:
        surf_dict[syms[i]].append(i)

    for i in bulk_indices:
        bulk_dict[syms[i]].append(i)

    # Get a list of bulk indices that CAN be swapped for each
    # unique surface element
    swap_dict = {}
    for surf_element in surf_dict:
        swap_dict[surf_element] = []
        for bulk_element in bulk_dict:
            if bulk_element == surf_element:
                continue
            swap_dict[surf_element] += bulk_dict[bulk_element]

    # Get a list of all swaps, taken based on their probability of happening        
    swap_list = []
    for i in range(min([len(surf_indices), len(bulk_indices)])):
        # First pick an element to swap based on the concentration in
        # the available surface swap atoms and atoms it can swap with
        surf_prob = np.asarray([len(surf_dict[element]) * len(swap_dict[element])
                                for element in surf_elements], dtype=float)
        if np.sum(surf_prob) == 0:
            break
        surf_prob /= np.sum(surf_prob)
        surf_element = np.random.choice(surf_elements, p=surf_prob)
        surf_index = np.random.choice(surf_dict[surf_element])

        # Pick a random bulk atom to swap to
        bulk_index = np.random.choice(swap_dict[surf_element])
        swap_list.append([surf_index, bulk_index])

        # Delete these indexes from the lists to avoid double swapping
        # either surface of bulk atoms
        surf_dict[surf_element].pop(surf_dict[surf_element].index(surf_index))
        for element in swap_dict:
            if bulk_index in swap_dict[element]:
                swap_dict[element].pop(swap_dict[element].index(bulk_index))

    if type(max_natoms) is float:
        max_natoms = int(max_natoms * len(swap_list))

    swap_list = swap_list[:max_natoms]

    # Finally swap the positions
    for surf_index, bulk_index in swap_list:
        surf_element = syms[surf_index]
        bulk_element = syms[bulk_index]
        syms[surf_index] = bulk_element
        syms[bulk_index] = surf_element

    individual.set_chemical_symbols(syms)

    return None
