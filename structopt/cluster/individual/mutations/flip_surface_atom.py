import random
import numpy as np

from structopt.tools import NeighborList
from structopt.tools import get_avg_radii

from ase.io import write

def flip_surface_atoms(individual, surf_CN=11, min_thickness=None):
    """Randomly moves atoms at the surface to other surface sites

    Parameters
    ----------
    individual : structopt.Individual object
        The individual object to be modified in place
    max_natoms : float or int
        if float, the maximum number of atoms that will be moved is 
        max_natoms*len(individual)
        if int, the maximum number of atoms that will be moved is max_natoms
        default: 0.20
    max_CN : int
        The coordination number cutoff for determining which atoms are surface atoms
        Any atoms with coordnation number at or above CN will not be considered as 
        surface.

    Output
    ------
    out : None
        Modifies individual in-place
    """

    if len(individual) == 0:
        return

    if min_thickness is None:
        min_thickness = get_avg_radii(individual)
    
    # Analyze the individual
    NNs = NeighborList(individual)
    CNs = [len(NN) for NN in NNs]
    
    # Get all surface atoms
    positions = individual.get_positions()
    surf_indices_CNs = [[i, CN] for i, CN in enumerate(CNs)
                        if CN <= surf_CN and CN > 2]
    surf_indices_CNs.sort(key=lambda i: i[1])
    surf_indices = list(zip(*surf_indices_CNs))[0]
    surf_positions = np.array([positions[i] for i in surf_indices])

    # Pair each surface atom with a surface atom on the other side
    # of the particle
    surf_xys = np.array([surf_positions[:,:2]])
    surf_xy_vecs = surf_xys - np.transpose(surf_xys, [1, 0, 2])
    surf_xy_dists = np.linalg.norm(surf_xy_vecs, axis=2)

    surf_zs = np.array(surf_positions[:,-1])
    surf_z_dists = surf_zs - surf_zs.T

    flips = (surf_z_dists > min_thickness).astype(float)
    flips[flips == 0] = np.inf
    
    
    surf_xy_dists *= np.inf
    
    # Get the average bond length of the particle
    chemical_symbols = individual.get_chemical_symbols()
    unique_symbols = set(chemical_symbols)
    atomlist = [[symbol, chemical_symbols.count(symbol)] for symbol in unique_symbols]
    avg_bond_length = get_avg_radii(atomlist) * 2

    # Set positions of a fraction of the surface atoms
    if type(max_natoms) is float:
        max_natoms = int(max_natoms * len(move_indices))
    move_natoms = random.randint(0, max_natoms)
    move_indices = move_indices[:move_natoms]
    add_indices = np.random.choice(len(add_positions), len(move_indices), replace=False)
    for move_index, add_index in zip(move_indices, add_indices):
        positions[move_index] = add_positions[add_index]

    individual.set_positions(positions)

    return move_indices
