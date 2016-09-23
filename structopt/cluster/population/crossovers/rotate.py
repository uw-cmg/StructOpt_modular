import random
import numpy as np
from ase import Atoms

from structopt.common.individual import Individual
from structopt.tools import root, single_core, parallel
from structopt.tools import random_three_vector

@single_core
def rotate(individual1, individual2, center_at_atom=True, repair_composition=True):
    """Rotates the two individuals around their centers of mass,
    splits them in half at the xy-plane, then splices them together.
    Maintains number of atoms.

    Parameters
    ----------
    individual1 : Individual
        The first parent
    individual2 : Individual 
        The second parent
    conserve_composition : bool 
        If True, conserves composition. For crossovers that create children
        with more (less) atoms, atoms are taken from (added to) the surface
        of the particle. For incorrect atomic ratios, atomic symbols are
        randomly interchanged throughout the particle

    Returns:
        Individual: The first child
        Individual: The second child

    The children are returned without indicies.
    """

    # Preserve starting conditions of individual
    ind1c = individual1.copy()
    ind2c = individual2.copy()

    # Translate individuals so COP is at (0, 0, 0)
    cop1 = ind1c.get_positions().mean(axis=0)
    cop2 = ind2c.get_positions().mean(axis=0)

    if center_at_atom:
        pos1 = ind1c.get_positions()
        dists1 = np.linalg.norm(pos1 - cop1, axis=1)
        cop1 = pos1[np.argmin(dists1)]

        pos2 = ind2c.get_positions()
        dists2 = np.linalg.norm(pos2 - cop2, axis=1)
        cop2 = pos2[np.argmin(dists2)]

    ind1c.translate(-cop1)
    ind2c.translate(-cop2)

    # Pick a random rotation angle and vector
    rot_vec = random_three_vector()
    rot_angle = random.uniform(0, 180) * np.pi / 180.0

    # Rotate both individuals
    ind1c.rotate(rot_vec, a=rot_angle, center=(0, 0, 0))
    ind2c.rotate(rot_vec, a=rot_angle, center=(0, 0, 0))

    # Create the children
    child1 = ind1c[ind1c.positions[:,2] >= 0] + ind2c[ind2c.positions[:,2] < 0]
    child2 = ind2c[ind2c.positions[:,2] >= 0] + ind1c[ind1c.positions[:,2] < 0]

    # Repair the children
    syms1, syms2 = ind1c.get_chemical_symbols(), ind2c.get_chemical_symbols()
    if repair_composition and sorted(syms1) == sorted(syms2):
        if len(syms1

    child1.rotate(rot_vec, a=-rot_angle, center=(0, 0, 0))
    child1.center()
    child2.rotate(rot_vec, a=-rot_angle, center=(0, 0, 0))
    child2.center()

    return child1, child2
