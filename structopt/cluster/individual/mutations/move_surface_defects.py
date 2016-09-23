import random
import numpy as np

from structopt.tools import CoordinationNumbers
from structopt.tools import get_avg_radii

from ase.io import write

def move_surface_defects(individual, surf_CN=11):
    """Moves atoms around on the surface based on coordination number
    Moves a surface atom with a low CN to an atom with a high CN

    Parameters
    ----------
    individual : structopt.Individual object
        The individual object to be modified in place
    surf_CN : int
        The maximum coordination number to considered a surface atom

    Output
    ------
    out : None
        Modifies individual in-place
    """

    if len(individual) == 0:
        return

    # Analyze the individual
    CNs = CoordinationNumbers(individual)
    
    # Get indices, CNs, and positions of all surface sites
    surf_indices_CNs = [[i, CN] for i, CN in enumerate(CNs) if CN <= surf_CN]
    if len(surf_indices_CNs) == 0:
        return 0
    surf_indices, surf_CNs = list(zip(*surf_indices_CNs))
    surf_indices = list(surf_indices)
    surf_CNs = list(surf_CNs)
    surf_positions = np.array([individual.positions[i] for i in surf_indices])

    # Get differences of CNs to find potentially good moves
    surf_CN_diffs = np.array([surf_CNs]) - np.array([surf_CNs]).T
    min_CN_diff = np.min(surf_CN_diffs) - 1
    surf_CN_diffs -= min_CN_diff
    unique_CN_diffs = np.array(list(set(surf_CN_diffs.flatten())))
    surf_CN_ps =  2.0 ** unique_CN_diffs
    surf_CN_ps /= np.sum(surf_CN_ps)
    surf_CN = np.random.choice(unique_CN_diffs, p=surf_CN_ps)

    switch_indices = np.argwhere(surf_CN_diffs == surf_CN)

    old_index, new_index = switch_indices[random.randint(0, len(switch_indices) - 1)]
    new_position = individual.positions[surf_indices[new_index]]
    
    # Get the average bond length of the particle
    chemical_symbols = individual.get_chemical_symbols()
    unique_symbols = set(chemical_symbols)
    atomlist = [[symbol, chemical_symbols.count(symbol)] for symbol in unique_symbols]
    avg_bond_length = get_avg_radii(atomlist) * 2

    # Choose sites as projections one bond length away from COP
    COP = surf_positions.mean(axis=0)
    vec = new_position - COP
    vec /= np.linalg.norm(vec)
    add_position = new_position + vec * avg_bond_length * 0.5

    individual[surf_indices[old_index]].position = add_position

    return 
