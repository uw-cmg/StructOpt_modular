import functools.wraps

import structopt.parameters
from . import LAMMPS


class Relaxations(object):
    """ """
    def __init__(self):
        self.parameters = structopt.parameters.relaxations
        self.relaxations = [getattr(self, name) for name in self.parameters.modules]

    @staticmethod
    @functools.wraps(LAMMPS)
    def LAMMPS(individual):
        return LAMMPS(individual)

