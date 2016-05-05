import random
import numpy as np


def rotate_atoms(individual, max_natoms=0.20):
    """Randomly rotates a number of random atoms within the individual (in place).
    
    Args:
        individual (Individual): an individual
        max_natoms (float or int): if float, the maximum number of atoms that will be rotated is max_natoms*len(individual)
                                   if int, the maximum number of atoms that will be rotated is max_natoms
                                   default: 0.20
    """
    if len(individual):
        if isinstance(max_natoms, float):
            max_natoms = int(len(individual)*max_natoms)
        max_natoms = max(max_natoms, 1)
        natoms_to_rotate = random.randint(1, max_natoms)
        atom_indices = list(range(len(individual)))
        random.shuffle(atom_indices)  # Using random.shuffle on the indices guarantees no duplicates
        atom_indices = atom_indices[:natoms_to_rotate]
        atoms = individual[atom_indices]  # This creates a copy of the atoms I think
        del individual[atom_indices]

        axis = random.choice(['x', '-x', 'y', '-y', 'z', '-z'])
        angle = random.uniform(30, 180)
        atoms.rotate(axis, a=angle, center='COM', rotate_cell=False)
        individual.extend(atoms)
    return None

