import functools
import random

import structopt
from .do_nothing import do_nothing
from .cost import cost


class Selections(object):
    """ """
    def __init__(self):
        self.parameters = structopt.parameters.selections
        self.selections = [getattr(self, name) for name in self.parameters.options]
        self.selected_selection = None

    def select_selection(self):
        self.selected_selection = random.choice(self.selections)

    def select(self, population, fits, nkeep):
        return self.selected_selection(population, fits, nkeep)

    @staticmethod
    @functools.wraps(do_nothing)
    def do_nothing(population, fits, nkeep):
        return do_nothing(population, fits, nkeep)

    @staticmethod
    @functools.wraps(cost)
    def cost(population, fits, nkeep):
        return cost(population, fits, nkeep)
