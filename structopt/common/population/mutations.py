import functools
import random

import structopt
from structopt.tools import root, single_core, parallel


class Mutations(object):
    """ """

    @single_core
    def __init__(self, parameters):
        self.parameters = parameters


    @single_core
    def mutate(self, population):
        for individual in population:
            individual.mutations.select_mutation()
            if individual.mutations.selected_mutation is not None:
                individual.mutate()
        return population


    @single_core
    def post_processing(self):
        pass

