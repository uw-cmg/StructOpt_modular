import logging
import functools
import importlib

import structopt
from . import LAMMPS
from . import hard_sphere_cutoff
from structopt.tools import root, single_core, parallel


class Relaxations(object):
    """Holds the parameters for each relaxation module and defines a utility function to run the relaxations for each relaxation module."""

    @single_core
    def __init__(self):
        self.parameters = structopt.parameters.relaxations
        self.modules = [globals()[module] for module in self.parameters.modules]


    @parallel
    def relax(self, population):
        """Relax the entire population using all the input relaxation methods.

        Args:
            population (Population): the population to relax
        """
        logger = logging.getLogger("default")
        to_relax = [individual for individual in population if individual._modified]
        logger.info("Relaxing individuals: {}".format(to_relax))
        for i, module in enumerate(self.modules):
            module.relax(population)

        return None


    @single_core
    def post_processing(self):
        pass

