import functools
import logging
import importlib
import numpy as np

import structopt
from . import LAMMPS, FEMSIM
from structopt.tools import root, single_core, parallel


class Fitnesses(object):
    """Holds the parameters for each fitness module and defines a utility function to compute the fitnesses for each fitness module."""

    @single_core
    def __init__(self):
        self.parameters = structopt.parameters.fitnesses

        self.modules = [globals()[module] for module in self.parameters.modules]


    @parallel
    def fitness(self, population):
        """Perform the fitness calculations on an entire population.

        Args:
            population (Population): the population to evaluate
        """
        fitnesses = np.zeros((len(population),), dtype=np.float)
        # Run each fitness module on the population
        for i, module in enumerate(self.modules):
            if structopt.parameters.globals.rank == 0:
                print("Running fitness {} on the entire population".format(module.__name__.split('.')[-1]))
            fits = module.fitness(population)

            # Calculate the full objective function with weights
            fits = np.multiply(fits, self.parameters.weights[i])
            fitnesses = np.add(fitnesses, fits)

        self.post_processing(fitnesses)
        return fitnesses


    @single_core
    def post_processing(self, fitnesses):
        logger = logging.getLogger("output")
        logger.info("Total fitnesses for the population: {} (rank {})".format(fitnesses, structopt.parameters.globals.rank))

