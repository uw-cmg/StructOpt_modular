import random
import numpy as np
from ase import Atom


def rotate_cluster(individual, max_natoms=0.20):
    """Randomly rotates a random cluster of atoms within the individual (in place).
    
    Args:
        individual (Individual): an individual
        max_natoms (float or int): if float, the maximum number of atoms that will be rotated is max_natoms*len(individual)
                                   if int, the maximum number of atoms that will be rotated is max_natoms
                                   default: 0.20
    """
    if len(individual):
        cell_max = np.maximum.reduce(individual.get_cell())
        cell_min = np.minimum.reduce(individual.get_cell())

        if isinstance(max_natoms, float):
            max_natoms = int(len(individual)*max_natoms)
        natoms_to_rotate = random.randint(1, max_natoms)

        point = (random.uniform(cell_min, cell_max),
                 random.uniform(cell_min, cell_max),
                 random.uniform(cell_min, cell_max))
        atom = Atom('Si', point)
        nearest_indices = individual.get_nearest_atom_indices(atom=atom, count=natoms_to_rotate)


        atoms = individual[nearest_indices]  # This creates a copy of the atoms I think
        del individual[nearest_indices]

        axis = random.choice(['x', '-x', 'y', '-y', 'z', '-z'])
        angle = random.uniform(30, 180)
        atoms.rotate(axis, a=angle, center='COM', rotate_cell=False)
        individual.extend(atoms)
    return None

