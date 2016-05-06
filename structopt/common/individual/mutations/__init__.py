import functools
import random
from itertools import accumulate
from bisect import bisect
from collections import defaultdict

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
        self.kwargs = defaultdict(dict)
        self.kwargs.update( {getattr(self, name): kwords for name, kwords in self.parameters.kwargs.items()} )

        self.mutations = {getattr(self, name): prob for name, prob in self.parameters.options.items()}
        total_probability = sum(self.mutations.values())
        self.mutations[None] = 1.0 - total_probability

        self.selected_mutation = None

    def select_mutation(self):
        # Implementation from https://docs.python.org/3/library/random.html -- Ctrl+F "weights"
        choices, weights = zip(*self.mutations.items())
        cumdist = list(accumulate(weights))
        x = random.random() * cumdist[-1]
        self.selected_mutation = choices[bisect(cumdist, x)]

    def mutate(self, individual):
        if self.selected_mutation is None:
            return individual
        else:
            individual._modified = True
            return self.selected_mutation(individual, **self.kwargs[self.selected_mutation])

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

