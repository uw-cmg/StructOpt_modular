import functools
import random

import structopt
from .do_nothing import do_nothing


class Predators(object):
    """ """
    def __init__(self):
        self.parameters = structopt.parameters.predators
        self.predators = [getattr(self, name) for name in self.parameters.options]
        self.selected_predator = None

    def select_predator(self):
        self.selected_predator = random.choice(self.predators)

    def kill(self, individual):
        return self.selected_predator(individual)

    @staticmethod
    @functools.wraps(do_nothing)
    def do_nothing(individual):
        return do_nothing(individual)

