import functools

import structopt
from . import LAMMPS, FEMSIM


class Fitnesses(object):
    """ """
    def __init__(self):
        self.parameters = structopt.parameters.fitnesses
        self.fitnesses = [getattr(self, name) for name in self.parameters.modules]

    def fitness(self, individual):
        return 0.0

    @staticmethod
    @functools.wraps(LAMMPS)
    def LAMMPS(individual):
        return LAMMPS(individual)

    @staticmethod
    @functools.wraps(FEMSIM)
    def FEMSIM(individual):
        return FEMSIM(individual)

