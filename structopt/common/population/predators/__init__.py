import functools
import random
from itertools import accumulate
from bisect import bisect

import structopt


class Predators(object):
    """ """
    def __init__(self):
        self.parameters = structopt.parameters.predators
        self.predators = {getattr(self, name): prob for name, prob in self.parameters.options.items()}
        total_probability = sum(self.predators.values())
        self.predators[None] = 1.0 - total_probability
        self.selected_predator = None

    def select_predator(self):
        # Implementation from https://docs.python.org/3/library/random.html -- Ctrl+F "weights"
        choices, weights = zip(*self.predators.items())
        cumdist = list(accumulate(weights))
        x = random.random() * cumdist[-1]
        self.selected_predator = choices[bisect(cumdist, x)]

    def kill(self, individual1, individual2):
        if self.selected_predator is None:
            return None, None
        else:
            return self.selected_predator(individual1, individual2)

