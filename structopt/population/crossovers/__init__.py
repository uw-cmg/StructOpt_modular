import functools
import random

import structopt
from .do_nothing import do_nothing


class Crossovers(object):
    """ """
    def __init__(self):
        self.parameters = structopt.parameters.crossovers
        self.crossovers = [getattr(self, name) for name in self.parameters.options]
        self.selected_crossover = None

    def select_crossover(self):
        self.selected_crossover = random.choice(self.crossovers)

    def crossover(self, individual):
        return self.selected_crossover(individual)

    @staticmethod
    @functools.wraps(do_nothing)
    def do_nothing(individual):
        return do_nothing(individual)

