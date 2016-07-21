import functools
import random
from itertools import accumulate
from bisect import bisect

import structopt
from structopt.tools import root, single_core, parallel
from .rotate import rotate


class Crossovers(object):
    """ """

    @single_core
    def __init__(self, parameters):
        self.parameters = parameters
        self.crossovers = {getattr(self, name): self.parameters[name]['probability'] for name in self.parameters}
        self.kwargs = {name: self.parameters[name]['kwargs'] for name in self.parameters}
        total_probability = sum(self.crossovers.values())
        self.crossovers[None] = 1.0 - total_probability
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
                kwargs = self.kwargs[self.selected_crossover.__name__.split('.')[-1]]
                child1, child2 = self.selected_crossover(individual1, individual2, **kwargs)
                if child1 is not None:
                    children.append(child1)
                if child2 is not None:
                    children.append(child2)
        return children


    @single_core
    def post_processing(self):
        pass

    @staticmethod
    @functools.wraps(rotate)
    def rotate(individual1, individual2):
        return rotate(individual1, individual2)

