import logging
import functools
import random
from itertools import accumulate
from bisect import bisect
from collections import defaultdict

import structopt
from structopt.tools import root, single_core, parallel

from .swap_positions import swap_positions
from .swap_species import swap_species
from .move_atoms import move_atoms
from .rotate_atoms import rotate_atoms
from .rotate_cluster import rotate_cluster
from .rotate_all import rotate_all
from .move_surface_atoms import move_surface_atoms

swap_positions.tag = 'SwPo'
swap_species.tag = 'SwSp'
move_atoms.tag = 'MoAt'
rotate_atoms.tag = 'RoAt'
rotate_cluster.tag = 'RoCl'
rotate_all.tag = 'RoAl'
move_surface_atoms.tag = 'MoSuAt'

class Mutations(object):
    """ """

    @single_core
    def __init__(self, parameters):
        # These variables never change
        self.parameters = parameters

        # self.mutations is a dictionary containing {function: probability} pairs
        self.mutations = {getattr(self, name): self.parameters[name]['probability'] for name in self.parameters}

        #self.kwargs is a dictionary containing {function: kwargs} pairs
        self.kwargs = {getattr(self, name): self.parameters[name]['kwargs'] for name in self.parameters}

        total_probability = sum(self.mutations.values())
        self.mutations[None] = 1.0 - total_probability

        # This parameter does not need to exist between generations
        self.selected_mutation = None


    @single_core
    def select_mutation(self):
        # Implementation from https://docs.python.org/3/library/random.html -- Ctrl+F "weights"
        choices, weights = zip(*self.mutations.items())
        cumdist = list(accumulate(weights))
        x = random.random() * cumdist[-1]
        self.selected_mutation = choices[bisect(cumdist, x)]


    @single_core
    def mutate(self, individual):
        logger = logging.getLogger("default")
        logger.info("Performing mutation {} on individual {}".format(self.selected_mutation, individual.id))
        print("Performing mutation {} on individual {}".format(self.selected_mutation, individual.id))
        if self.selected_mutation is None:
            return individual
        else:
            individual._relaxed = False
            individual._fitted = False
            kwargs = self.kwargs[self.selected_mutation]
            ret = self.selected_mutation(individual, **kwargs)
            self.post_processing(individual)
            return ret


    @single_core
    def post_processing(self, individual):
        individual.mutation_tag = 'm{tag}'.format(tag=self.selected_mutation.tag)


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


    @staticmethod
    @functools.wraps(rotate_all)
    def rotate_all(individual, vector=None, angle=None, center=None):
        return rotate_all(individual, vector, angle, center)

    @staticmethod
    @functools.wraps(move_surface_atoms)
    def move_surface_atoms(individual, max_natoms=0.2, move_CN=9, surf_CN=11):
        return move_surface_atoms(individual, max_natoms, move_CN, surf_CN)
    

