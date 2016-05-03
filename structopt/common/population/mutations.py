import functools
import random

import structopt


class Mutations(object):
    """ """
    def __init__(self):
        self.parameters = structopt.parameters.mutations
        self.mutations = self.parameters.options  # [getattr(self, name) for name in self.parameters.options] TODO: this is not consistent between the population and individual level. how do I handle this?
        self.selected_mutation = None

    def select_mutation(self):
        self.selected_mutation = random.choice(self.mutations)

    def mutate(self, population):
        for individual in population:
            individual.mutations.select_mutation(self.selected_mutation)
            individual.mutate()
        return population

