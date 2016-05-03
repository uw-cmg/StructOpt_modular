import functools
import random

import structopt
from .add_atoms import add_atoms
from .remove_atoms import remove_atoms
from .remove_surface_atoms import remove_surface_atoms
from .lattice_alteration import lattice_alteration
from .lattice_alteration_group import lattice_alteration_group
from .rotation import rotation
from .rotation_geo import rotation_geo


class Mutations(object):
    """ """
    def __init__(self):
        self.parameters = structopt.parameters.mutations
        self.mutations = [getattr(self, name) for name in self.parameters.options]
        self.selected_mutation = None

    def select_mutation(self, name=None):
        if name is not None:
            self.selected_mutation = getattr(self, name)
        else:
            self.selected_mutation = random.choice(self.mutations)

    def mutate(self, individual):
        return self.selected_mutation(individual)

    @staticmethod
    @functools.wraps(add_atoms)
    def add_atoms(individual):
        return add_atoms(individual)

    @staticmethod
    @functools.wraps(remove_atoms)
    def remove_atoms(individual):
        return remove_atoms(individual)

    @staticmethod
    @functools.wraps(remove_surface_atoms)
    def remove_surface_atoms(individual):
        return remove_surface_atoms(individual)

    @staticmethod
    @functools.wraps(lattice_alteration)
    def lattice_alteration(individual):
        return lattice_alteration(individual)

    @staticmethod
    @functools.wraps(lattice_alteration_group)
    def lattice_alteration_group(individual):
        return lattice_alteration_group(individual)

    @staticmethod
    @functools.wraps(rotation)
    def rotation(individual):
        return rotation(individual)

    @staticmethod
    @functools.wraps(rotation_geo)
    def rotation_geo(individual):
        return rotation_geo(individual)

