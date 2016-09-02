import functools
import structopt.common.individual.mutations

from .move_atoms import move_atoms
from .move_surface_atoms import move_surface_atoms
from .move_atoms_group import move_atoms_group
from .rotate_cluster import rotate_cluster
from .twist import twist
from .swap_core_shell import swap_core_shell
from .rich2poor import rich2poor
from .poor2rich import poor2rich

move_surface_atoms.tag = 'MoSuAt'
move_atoms.tag = 'MoAt'
move_atoms_group.tag = 'MoAtGr'
rotate_cluster.tag = 'RoCl'
twist.tag = 'Twist'
swap_core_shell.tag = 'SwCoSh'
rich2poor.tag = 'Ri2Po'

class Mutations(structopt.common.individual.mutations.Mutations):

    @staticmethod
    @functools.wraps(move_surface_atoms)
    def move_surface_atoms(individual, max_natoms=0.2, move_CN=9, surf_CN=11):
        return move_surface_atoms(individual, max_natoms, move_CN, surf_CN)
    
    @staticmethod
    @functools.wraps(move_atoms)
    def move_atoms(individual, max_natoms=0.20):
        return move_atoms(individual, max_natoms)

    @staticmethod
    @functools.wraps(move_atoms_group)
    def move_atoms_group(individual, max_natoms=0.20):
        return move_atoms_group(individual, max_natoms)

    @staticmethod
    @functools.wraps(rotate_cluster)
    def rotate_cluster(individual, max_natoms=0.20):
        return rotate_cluster(individual, max_natoms)

    @staticmethod
    @functools.wraps(twist)
    def twist(individual, max_radius=0.90):
        return twist(individual, max_radius)
    
    @staticmethod
    @functools.wraps(swap_core_shell)
    def swap_core_shell(individual, max_natoms=0.2, surf_CN=11):
        return twist(individual, max_natoms, surf_CN)
