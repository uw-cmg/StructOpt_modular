import functools
import random

import structopt


class Mutations(object):
    """ """
    def __init__(self):
        self.parameters = structopt.parameters.mutations

    def mutate(self, population):
        for individual in population:
            individual.mutations.select_mutation()
            if individual.selected_mutation is not None:
                individual.mutate()
        return population

