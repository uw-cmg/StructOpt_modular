import logging
import subprocess

import structopt
from structopt.tools import root, single_core, parallel
from structopt.tools.parallel import allgather
import numpy as np


@parallel
def fitness(population, parameters):
    """Perform the LAMMPS fitness calculation on an entire population.

    Args:
        population (Population): the population to evaluate
    """
    if parameters.use_mpi4py:
        logger = logging.getLogger('by-rank')
    else:
        logger = logging.getLogger('output')

    to_fit = [individual for individual in population if not individual._fitted]
    if parameters.use_mpi4py:
        ncores = logging.parameters.ncores
    else:
        ncores = 1
    rank = logging.parameters.rank

    individuals_per_core = {r: [] for r in range(ncores)}
    for i, individual in enumerate(to_fit):
        individuals_per_core[i % ncores].append(individual)

    for individual in individuals_per_core[rank]:
        print("Running LAMMPS fitness evaluation on individual {}".format(individual.index))
        energy = individual.fitnesses.LAMMPS.fitness(individual, population.generation)
        individual.LAMMPS = energy
        logger.info('Individual {0} after LAMPPS evaluation has energy {1}'.format(individual.index, energy))

    fits = [individual.LAMMPS for individual in population]
    positions_per_core = {rank: [population.position(individual) for individual in individuals] for rank, individuals in individuals_per_core.items()}
    if parameters.use_mpi4py:
        fits = allgather(fits, positions_per_core)

    # Save the fitness value for the module to each individual after they have been allgathered
    for i, fit in enumerate(fits):
        population.get_by_position(i).LAMMPS = fit

    return [individual.LAMMPS for individual in population]

