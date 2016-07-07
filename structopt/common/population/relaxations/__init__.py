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
    def __init__(self, parameters=None):
        self.parameters = parameters or structopt.parameters.relaxations
        self.modules = [globals()[module] for module in self.parameters.modules]


    @parallel
    def relax(self, population):
        """Relax the entire population using all the input relaxation methods.

        Args:
            population (Population): the population to relax
        """
        logger = logging.getLogger("default")
        to_relax = [individual for individual in population if not individual._relaxed]
        logger.info("Found {} individuals to relax on core {}: {}".format(len(to_relax), structopt.parameters.globals.rank, to_relax))
        #print("Found {} individuals to relax on core {}: {}".format(len(to_relax), structopt.parameters.globals.rank, to_relax))
        for i, module in enumerate(self.modules):
            if structopt.parameters.globals.rank == 0:
                print("Running relaxation {} on the entire population".format(module.__name__.split('.')[-1]))
            module.relax(population)

        for individual in population:
            individual._relaxed = True
        return None


    @single_core
    def post_processing(self):
        pass

