import functools
import random

import structopt
from .do_nothing import do_nothing


class Selections(object):
    """ """
    def __init__(self):
        self.parameters = structopt.parameters.selections
        self.selections = [getattr(self, name) for name in self.parameters.options]
        self.selected_selection = None

    def select_selection(self):
        self.selected_selection = random.choice(self.selections)

    def select(self, individual):
        return self.selected_selection(individual)

    @staticmethod
    @functools.wraps(do_nothing)
    def do_nothing(individual):
        return do_nothing(individual)

