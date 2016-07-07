import functools
import random
from itertools import accumulate
from bisect import bisect

import structopt
from .random_selection import random_selection
from .rank import rank
from structopt.tools import root, single_core, parallel


class Selections(object):
    """ """

    @single_core
    def __init__(self, parameters=None):
        self.parameters = parameters or structopt.parameters.selections
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
    def select(self, population, fits):
        pairs = self.selected_selection(population=population,
                                        fits=fits,
                                        prob=self.parameters.crossover_probability)
        self.post_processing(pairs)
        return pairs


    @single_core
    def post_processing(self, pairs):
        pass


    @staticmethod
    @functools.wraps(random_selection)
    def random_selection(population, fits, prob):
        return random_selection(population, fits, prob=prob)

    @staticmethod
    @functools.wraps(rank)
    def rank(population, fits, prob):
        return rank(population, fits, prob=prob)
