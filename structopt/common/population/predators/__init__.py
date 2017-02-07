import functools
import random
from itertools import accumulate
from bisect import bisect
import numpy as np

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
    def kill(self, population, nkeep, keep_best=True):
        """Removes some individuals from the population.

        Parameters
        ----------
        nkeep : int
            The number of individuals to keep.
        keep_best : bool
            If set to True, the best individual is always included in the
            following generation.

        Returns
        -------
            The individuals that were removed from the population.
        """
        if self.selected_predator is None:
            return
        assert nkeep <= len(population)
        if nkeep <= len(population)
            return

        fits = {individual.id: individual.fitness for individual in population}
        if keep_best:
            best_id = min(fits, key=fits.get)
            fits.pop(best_id)
            nkeep -= 1

        kwargs = self.kwargs[self.selected_predator]
        to_keep = self.selected_predator(fits=fits, nkeep=nkeep, **kwargs)
        if keep_best:
            try:
                to_keep.append(best_id)
            except AttributeError:
                to_keep = np.append(to_keep, best_id)

        killed = [individual for individual in population if individual.id not in to_keep]
        new_population = [population[id] for id in to_keep]
        population.replace(new_population)

        return killed

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

