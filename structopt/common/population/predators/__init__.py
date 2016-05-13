import functools
import random
from itertools import accumulate
from bisect import bisect

import structopt
from structopt.tools import root, single_core, parallel


class Predators(object):
    """ """

    @single_core
    def __init__(self):
        self.parameters = structopt.parameters.predators
        self.predators = {getattr(self, name): prob for name, prob in self.parameters.options.items()}
        total_probability = sum(self.predators.values())
        self.predators[None] = 1.0 - total_probability
        self.selected_predator = None


    @single_core
    def select_predator(self):
        # Implementation from https://docs.python.org/3/library/random.html -- Ctrl+F "weights"
        choices, weights = zip(*self.predators.items())
        cumdist = list(accumulate(weights))
        x = random.random() * cumdist[-1]
        self.selected_predator = choices[bisect(cumdist, x)]


    @single_core
    def kill(self, population, fits):
        pass # TODO


    @single_core
    def post_processing(self):
        pass

