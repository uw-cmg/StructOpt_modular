import functools
import structopt.common.individual.mutations

from .move_atoms import move_atoms
from .move_surface_atoms import move_surface_atoms
from .move_atoms_group import move_atoms_group

move_surface_atoms.tag = 'MoSuAt'
move_atoms.tag = 'MoAt'
move_atoms_group.tag = 'MoAtGr'

class Mutations(structopt.common.individual.mutations.Mutations):

    @staticmethod
    @functools.wraps(move_surface_atoms)
    def move_surface_atoms(individual, max_natoms=0.2, move_CN=9, surf_CN=11):
        return move_surface_atoms(individual, max_natoms, move_CN, surf_CN)
    
    @staticmethod
    @functools.wraps(move_atoms)
    def move_atoms(individual, max_natoms=0.20):
        return move_atoms(individual, max_natoms=0.20)

    @staticmethod
    @functools.wraps(move_atoms_group)
    def move_atoms_group(individual, max_natoms=0.20):
        return move_atoms_group(individual, max_natoms=0.20)
