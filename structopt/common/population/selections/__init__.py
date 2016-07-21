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
    def __init__(self, parameters, crossover_probability):
        self.parameters = parameters
        self.selections = {getattr(self, name): self.parameters[name]['probability'] for name in self.parameters}
        self.kwargs = {name: self.parameters[name]['kwargs'] for name in self.parameters}
        assert sum(self.selections.values()) == 1.0
        self.selected_selection = None
        self.crossover_probability = crossover_probability


    @single_core
    def select_selection(self):
        # Implementation from https://docs.python.org/3/library/random.html -- Ctrl+F "weights"
        choices, weights = zip(*self.selections.items())
        cumdist = list(accumulate(weights))
        x = random.random() * cumdist[-1]
        self.selected_selection = choices[bisect(cumdist, x)]


    @single_core
    def select(self, population, fits):
        kwargs = self.kwargs[self.selected_selection.__name__.split('.')[-1]]
        pairs = self.selected_selection(population=population,
                                        fits=fits,
                                        prob=self.crossover_probability,
                                        **kwargs)
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
