import functools
import random
from itertools import accumulate
from bisect import bisect

import structopt
from structopt.tools import root, single_core, parallel

from .rotate import rotate

rotate.tag = 'Ro'


class Crossovers(object):
    """ """

    @single_core
    def __init__(self, parameters):
        self.parameters = parameters

        # self.crossovers is a dictionary containing {function: probability} pairs
        self.crossovers = {getattr(self, name): self.parameters[name]['probability'] for name in self.parameters}

        #self.kwargs is a dictionary containing {function: kwargs} pairs
        self.kwargs = {getattr(self, name): self.parameters[name]['kwargs'] for name in self.parameters}

        self.total_probability = sum(self.crossovers.values())
        self.crossovers[None] = 1.0 - self.total_probability
        self.selected_crossover = None


    @single_core
    def select_crossover(self):
        # Implementation from https://docs.python.org/3/library/random.html -- Ctrl+F "weights"
        choices, weights = zip(*self.crossovers.items())
        cumdist = list(accumulate(weights))
        x = random.random() * cumdist[-1]
        self.selected_crossover = choices[bisect(cumdist, x)]


    @single_core
    def crossover(self, pairs):
        children = []
        for individual1, individual2 in pairs:
            self.select_crossover()  # Choose a new crossover to perform for every pair
            if self.selected_crossover is not None:
                print("Performing crossover {} on individuals {} and {}".format(self.selected_crossover, individual1, individual2))
                kwargs = self.kwargs[self.selected_crossover]
                child1, child2 = self.selected_crossover(individual1, individual2, **kwargs)
                if child1 is not None:
                    children.append(child1)
                if child2 is not None:
                    children.append(child2)

                self.post_processing((individual1, individual2), (child1, child2))
        return children


    @single_core
    def post_processing(self, parent_pair, child_pair):
        parent1, parent2 = parent_pair
        child1, child2 = child_pair
        for child in child_pair:
            if child is not None:
                child.crossover_tag = '({parent1}+{parent2}){tag}'.format(parent1=parent1.index, parent2=parent2.index, tag=self.selected_crossover.tag)


    @staticmethod
    @functools.wraps(rotate)
    def rotate(individual1, individual2, conserve_composition=True):
        return rotate(individual1, individual2, conserve_composition)

