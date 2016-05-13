import functools
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


    @single_core
    def fitness(self, population):
        fitnesses = np.zeros((len(population),), dtype=np.float)
        for i, module in enumerate(self.modules):
            fits = module.fitness(population)

            # Save the fitness value for the module to each individual
            for j, individual in enumerate(population):
                setattr(individual, self.parameters.modules[i], fits[j])

            fits = np.multiply(fits, self.parameters.weights[i])
            fitnesses = np.add(fitnesses, fits)

        # Set each individual to unmodified so that the fitnesses don't need to be recalculated
        for individual in population:
            individual._modifed = False

        return fitnesses


    @single_core
    def post_processing(self):
        pass

