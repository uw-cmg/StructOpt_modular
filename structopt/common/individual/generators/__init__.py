import functools
import numpy as np
from importlib import import_module

import structopt
from .read_xyz import read_xyz
from structopt.tools import root, single_core, parallel

@single_core
def generate(individual, **kwargs):
    """ Uses the relevant parameters from structopt to intialize the 
    input Individual by modifying it in-place.

        Args:
            individual (Individual): an Individual that is uninitialized
            **kwargs: keyword arguments for either ase.Atoms 
                      or a different generator function
    """

    if not kwargs:
        return None

    # If we are reading from a file, load the atoms and return None
    if 'filenames' in kwargs or 'filename' in kwargs:
        if 'filenames' in kwargs:
            filename = kwargs['filenames'][individual.index]
        else:
            filename = kwargs['filename']
        atoms = read_xyz(filename)
        individual.extend(atoms)

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

    # Else we are calling a function. The name of the function is
    # the kwarg. It's parameters are given as another dictionary object

    # Currently we only support one type of generator
    # Remove this once multiple generators are supported
    assert len(kwargs) == 1

    # Get the generator object
    name = list(kwargs.keys())[0]
    module = 'structopt.common.individual.generators.{}'.format(name.lower())
    module = import_module(module)
    Generator = getattr(module, name.title())

    # Get an atoms object from the generator and modify the individual
    atomlist = structopt.parameters.generators.atomlist
    params = kwargs[name]
    generator = Generator(atomlist, **params)
    atoms = generator.generate()
    individual.extend(atoms)

    # Most of the time the generated atoms are centerd at the origin
    # Center them in the middle of the cell. LAMMPS doesn't wrap around
    individual.set_cell(structopt.parameters.generators.cell)
    center = np.sum(atoms.get_cell(), axis=0) * 0.5
    individual.translate(center)

    return None

