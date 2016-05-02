import functools.wraps

import structopt.parameters
from . import add_atom, remove_atom, remove_surface_atom


class Mutations(object):
    """ """
    def __init__(self):
        self.parameters = structopt.parameters.mutations
        self.mutations = [getattr(self, name) for name in self.parameters.options]
        self.selected_mutation = None

    def select_mutation(self):
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

