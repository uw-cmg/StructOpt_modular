import numpy as np

from structopt.tools import NeighborList
from structopt.tools import get_avg_radii

def move_surface_atoms(individual, max_natoms=0.2, move_CN=9, surf_CN=11):
    """Randomly moves atoms at the surface to other surface sites

    Parameters
    ----------
    individual : structopt.Individual object
        The individual object to be modified in place
    max_natoms : float or int
        if float, the maximum number of atoms that will be moved is max_natoms*len(individual)
        if int, the maximum number of atoms that will be moved is max_natoms
        default: 0.20
    max_CN : int
        The coordination number cutoff for determining which atoms are surface atoms
        Any atoms with coordnation number at or above CN will not be considered as surface.

    Output
    ------
    out : None
        Modifies individual in-place
    """

    if not len(individual):
        return

    # Analyze the individual
    NNs = NeighborList(individual)
    CNs = [len(NN) for NN in NNs]
    
    # Get indexes of atoms considered to be moved
    move_indexes_CNs = [[i, CN] for i, CN in enumerate(CNs) if CN < move_CN]
    move_indexes_CNs.sort(key=lambda i: i[1])
    move_indexes = list(zip(*move_indexes_CNs))[0]

    # Get surface sites to move atoms to
    # First get all surface atoms
    positions = individual.get_positions()
    surf_indexes_CNs = [[i, CN] for i, CN in enumerate(CNs) if CN < surf_CN]
    surf_indexes_CNs.sort(key=lambda i: i[1])
    surf_indexes = list(zip(*surf_indexes_CNs))[0]
    surf_positions = np.array([positions[i] for i in surf_indexes])

    # Get the average bond length of the particle
    chemical_symbols = individual.get_chemical_symbols()
    unique_symbols = set(chemical_symbols)
    atomlist = [[symbol, chemical_symbols.count(symbol)] for symbol in unique_symbols]
    avg_bond_length = get_avg_radii(atomlist) * 2

    # Choose sites as projections one bond length away from COM
    COM = individual.get_center_of_mass()    
    vec = surf_positions - COM
    vec /= np.array([np.linalg.norm(vec, axis=1)]).T
    add_positions = surf_positions + vec * avg_bond_length

    # Set positions of a fraction of the surface atoms
    move_indexes = move_indexes[:int(max_natoms * len(move_indexes))]
    add_indexes = np.random.choice(len(add_positions), len(move_indexes), replace=False)
    for move_index, add_index in zip(move_indexes, add_indexes):
        positions[move_index] = add_positions[add_index]

    individual.set_positions(positions)

    return None

