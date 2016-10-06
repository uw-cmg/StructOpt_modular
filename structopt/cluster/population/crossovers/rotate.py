import random
import numpy as np
from ase import Atoms

from structopt.common.individual import Individual
from structopt.tools import root, single_core, parallel
from structopt.tools import random_three_vector
from structopt.tools import repair_cluster

import time

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

    # Translate individuals so COP is at (0, 0, 0)
    cop1 = individual1.get_positions().mean(axis=0)
    cop2 = individual2.get_positions().mean(axis=0)

    if center_at_atom:
        pos1 = individual1.get_positions()
        dists1 = np.linalg.norm(pos1 - cop1, axis=1)
        cop1 = pos1[np.argmin(dists1)]

        pos2 = individual2.get_positions()
        dists2 = np.linalg.norm(pos2 - cop2, axis=1)
        cop2 = pos2[np.argmin(dists2)]

    individual1.translate(-cop1)
    individual2.translate(-cop2)

    # Pick a random rotation angle and vector
    rot_vec = random_three_vector()
    rot_angle = random.uniform(0, 180) * np.pi / 180.0

    # Rotate both individuals
    individual1.rotate(rot_vec, a=rot_angle, center=(0, 0, 0))
    individual2.rotate(rot_vec, a=rot_angle, center=(0, 0, 0))

    # Create the children
    child1 = Atoms()
    child2 = Atoms()
    child1.extend(individual1[individual1.positions[:,2] >= 0])
    child1.extend(individual2[individual2.positions[:,2] < 0])
    child2.extend(individual2[individual2.positions[:,2] >= 0])
    child2.extend(individual1[individual1.positions[:,2] < 0])

    # Repair the children
    if repair_composition:
        syms = individual1.get_chemical_symbols()
        atomlist = [[sym, syms.count(sym)] for sym in set(syms)]
        repair_cluster(child1, atomlist)
        repair_cluster(child2, atomlist)

    # Reorient the children
    child1.rotate(rot_vec, a=-rot_angle, center=(0, 0, 0))
    child1.center()
    child2.rotate(rot_vec, a=-rot_angle, center=(0, 0, 0))
    child2.center()

    # Reorient the parents
    individual1.rotate(rot_vec, a=-rot_angle, center=(0, 0, 0))
    individual1.translate(cop1)
    individual2.rotate(rot_vec, a=-rot_angle, center=(0, 0, 0))
    individual2.translate(cop2)

    t0 = time.time()
    full_child1 = individual1.copy(include_atoms=False)
    full_child1.extend(child1)
    full_child2 = individual2.copy(include_atoms=False)
    full_child2.extend(child2)
    print(time.time() - t0)
    
    return full_child1, full_child2
