import functools
import numpy as np

import structopt
from structopt.common.individual import Individual
from .read_xyz import read_xyz
from structopt.tools import root, single_core, parallel

@single_core
def generate(individual, **kwargs):
    """ Uses the relevant parameters from structopt to intialize the input Individual by modifying it in-place.

        Args:
            individual (Individual): an Individual that is uninitialized
            **kwargs: keyword arguments for either ase.Atoms or a different generator function
    """
    if 'filenames' in kwargs:
        filename = kwargs['filenames'][individual.index]
        atoms = read_xyz(filename)
        individual.extend(atoms)
        #individual.set_atomic_numbers(atoms.get_atomic_numbers())
        #individual.set_charges(atoms.get_charges())
        #individual.set_chemical_symbols(atoms.get_chemical_symbols())
        #individual.set_initial_magnetic_moments(atoms.get_initial_magnetic_moments())
        #individual.set_masses(atoms.get_masses())
        #individual.set_momenta(atoms.get_momenta())
        #individual.set_positions(atoms.get_positions())
        #individual.set_scaled_positions(atoms.get_scaled_positions())
        #individual.set_tags(atoms.get_tags())
        #individual.set_velocities(atoms.get_velocities())

        individual.set_pbc(True)
        with open(filename) as of:
            of.readline()  # number of atoms
            sizes = of.readline()  # comment == box size
            sizes = sizes.split()[:3]
            sizes = [float(x) for x in sizes]
            cell = np.identity(3)
            np.fill_diagonal(cell, sizes)
            individual.set_cell(cell)

    return None

