import functools
import importlib

import structopt
from . import LAMMPS
from . import hard_sphere_cutoff


class Relaxations(object):
    """Holds the parameters for each relaxation module and defines a utility function to run the relaxations for each relaxation module."""
    def __init__(self):
        self.parameters = structopt.parameters.relaxations
        self.modules = [globals()[module] for module in self.parameters.modules]


    def relax(self, population):
        for i, module in enumerate(self.modules):
            population = module.relax(population)

        return None

    def post_processing(self):
        pass

