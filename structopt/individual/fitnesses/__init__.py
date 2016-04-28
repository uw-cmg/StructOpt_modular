import functools.wraps

import structopt.parameters
from . import LAMMPS


class Fitnesses(object):
    """ """
    def __init__(self):
        self.parameters = structopt.parameters.fitnesses
        self.fitnesses = [getattr(self, name) for name in self.parameters.modules]

    @staticmethod
    @functools.wraps(LAMMPS)
    def LAMMPS(individual):
        return LAMMPS(individual)

