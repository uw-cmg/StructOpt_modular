import random
import numpy as np

from structopt.tools import NeighborList
from structopt.tools import get_avg_radii

from ase.io import write

def flip_surface_atom(individual, surf_CN=11, min_thickness=None):
    """Randomly "flips" an atom from one side of the particle to the other side
    of the particle, ideally making little changes to the STEM fitness

    Parameters
    ----------
    individual : structopt.Individual object
        The individual object to be modified in place
    surf_CN : int
        The maximum coordination number for determining an atom as a surface atom
    min_thickness : float
        The minimum distance to consider an atom as on the other side
        of the particle.

    Output
    ------
    out : None
        Modifies individual in-place
    """

    if len(individual) == 0:
        return

    if min_thickness is None:
        min_thickness = get_avg_radii(individual) * 2
    
    # Analyze the individual
    NNs = NeighborList(individual)
    CNs = [len(NN) for NN in NNs]
    
    # Get all surface atoms
    positions = individual.get_positions()
    surf_indices_CNs = [[i, CN] for i, CN in enumerate(CNs)
                        if CN <= surf_CN and CN > 2]
    surf_indices, surf_CNs = list(zip(*surf_indices_CNs))

    # Pair each surface atom with a surface atom on the other side
    # of the particle
    surf_positions = np.array([positions[i] for i in surf_indices])    
    surf_xys = np.array([surf_positions[:,:2]])
    surf_xy_vecs = surf_xys - np.transpose(surf_xys, [1, 0, 2])
    surf_xy_dists = np.linalg.norm(surf_xy_vecs, axis=2)

    surf_zs = np.array([surf_positions[:,-1]])
    surf_z_dists = np.absolute(surf_zs - surf_zs.T)

    flips = (surf_z_dists > min_thickness).astype(float)
    flips[flips == 0] = np.inf
    
    surf_xy_dists *= flips
    np.fill_diagonal(surf_xy_dists, np.inf)
    flips = [np.argmin(dists) for dists in surf_xy_dists]
    flip_indices = [[surf_indices[i], surf_indices[j]] for i, j in enumerate(flips)]
    flip_CNs_diffs = [CNs[i] - CNs[j] for i, j in flip_indices]

    # Simple rank based selection 

    

    individual.set_positions(positions)

    return
