import functools
import random
from itertools import accumulate
from bisect import bisect

from structopt.tools import root, single_core, parallel
from .best import best
from .roulette import roulette
from .tournament import tournament
from .rank import rank
from .fuss import fuss


class Predators(object):
    """ """

    @single_core
    def __init__(self, parameters):
        self.parameters = parameters
        self.predators = {getattr(self, name): self.parameters[name]['probability'] for name in self.parameters}
        self.kwargs = {getattr(self, name): self.parameters[name]['kwargs'] for name in self.parameters}
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
    def kill(self, population, fits, nkeep):
        if self.selected_predator is not None:
            kwargs = self.kwargs[self.selected_predator]
            self.selected_predator(population=population, fits=fits, nkeep=nkeep, **kwargs)

    @single_core
    def post_processing(self):
        pass

    @staticmethod
    @functools.wraps(best)
    def best(population, fits, nkeep):
        return best(population, fits, nkeep)

    @staticmethod
    @functools.wraps(roulette)
    def roulette(population, fits, nkeep, keep_best=True, T=None):
        return roulette(population, fits, nkeep, keep_best, T)

    @staticmethod
    @functools.wraps(tournament)
    def tournament(population, fits, nkeep, tournament_size=5, keep_best=True):
        return tournament(population, fits, nkeep, tournament_size, keep_best)

    @staticmethod
    @functools.wraps(rank)
    def rank(population, fits, nkeep, p_min=None):
        return rank(population, fits, nkeep, p_min)

    @staticmethod
    @functools.wraps(fuss)
    def fuss(population, fits, nkeep, nbest=1, fusslimit=10):
        return fuss(population, fits, nkeep, nbest=1, fusslimit=10)

