import functools
import random

import structopt
from .swap_positions import swap_positions
from .swap_species import swap_species
from .move_atoms import move_atoms
from .rotate_atoms import rotate_atoms
from .rotate_cluster import rotate_cluster


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
    @functools.wraps(swap_positions)
    def swap_positions(individual):
        return swap_positions(individual)

    @staticmethod
    @functools.wraps(swap_species)
    def swap_species(individual):
        return swap_species(individual)

    @staticmethod
    @functools.wraps(move_atoms)
    def move_atoms(individual):
        return move_atoms(individual)

    @staticmethod
    @functools.wraps(rotate_atoms)
    def rotate_atoms(individual):
        return rotate_atoms(individual)

    @staticmethod
    @functools.wraps(rotate_cluster)
    def rotate_cluster(individual):
        return rotate_cluster(individual)

