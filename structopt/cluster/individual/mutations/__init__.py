import functools
import structopt.common.individual.mutations

from .move_atoms import move_atoms
from .move_surface_atoms import move_surface_atoms
from .move_surface_STEM import move_surface_STEM
from .move_atoms_group import move_atoms_group
from .rotate_cluster import rotate_cluster
from .twist import twist
from .swap_core_shell import swap_core_shell
from .rich2poor import rich2poor
from .poor2rich import poor2rich
from .flip_surface_atom import flip_surface_atom
from .permute_column_STEM import permute_column_STEM
from .permute_column_surface import permute_column_surface

move_surface_atoms.tag = 'MoSuAt'
move_surface_STEM.tag = 'MoSuSTEM'
move_atoms.tag = 'MoAt'
move_atoms_group.tag = 'MoAtGr'
rotate_cluster.tag = 'RoCl'
twist.tag = 'Twist'
swap_core_shell.tag = 'SwCoSh'
rich2poor.tag = 'Ri2Po'
poor2rich.tag = 'Pr2Ri'
permute_column_surface.tag = 'PeCoSu'

class Mutations(structopt.common.individual.mutations.Mutations):

    @staticmethod
    @functools.wraps(move_surface_atoms)
    def move_surface_atoms(individual, max_natoms=0.2, move_CN=8, surf_CN=11):
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
        return swap_core_shell(individual, max_natoms, surf_CN)

    @staticmethod
    @functools.wraps(move_surface_STEM)
    def move_surface_STEM(individual, STEM_parameters, move_CN=11, surf_CN=11,
                          filter_size=1, move_cutoff=0.5, surf_cutoff=0.5,
                          max_cutoff=0.5, min_cutoff=0.5):
        return move_surface_STEM(individual, STEM_parameters, move_CN, surf_CN,
                                 filter_size, move_cutoff, surf_cutoff,
                                 max_cutoff, min_cutoff)

    @staticmethod
    @functools.wraps(rich2poor)
    def rich2poor(individual, max_natoms=0.05, surf_CN=11, factor=1.1):
        return rich2poor(individual, max_natoms, surf_CN, factor)

    @staticmethod
    @functools.wraps(permute_column_surface)
    def permute_column_surface(individual, STEM_parameters, filter_size=0.5,
                               column_cutoff=0.5):
        return return permute_column_surface(individual, STEM_parameters,
                                             filter_size, column_cutoff)
