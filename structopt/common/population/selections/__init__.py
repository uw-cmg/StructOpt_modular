import functools
import random
from itertools import accumulate
from bisect import bisect

import structopt
from .cost import cost
from structopt.tools import root, single_core, parallel


class Selections(object):
    """ """

    @single_core
    def __init__(self):
        self.parameters = structopt.parameters.selections
        self.selections = {getattr(self, name): prob for name, prob in self.parameters.options.items()}
        assert sum(self.selections.values()) == 1.0
        self.selected_selection = None


    @single_core
    def select_selection(self):
        # Implementation from https://docs.python.org/3/library/random.html -- Ctrl+F "weights"
        choices, weights = zip(*self.selections.items())
        cumdist = list(accumulate(weights))
        x = random.random() * cumdist[-1]
        self.selected_selection = choices[bisect(cumdist, x)]


    @single_core
    def select(self, population, fits, nkeep):
        return self.selected_selection(population, fits, nkeep)


    @single_core
    def post_processing(self):
        pass

    @staticmethod
    @functools.wraps(cost)
    def cost(population, fits, nkeep):
        return cost(population, fits, nkeep)

