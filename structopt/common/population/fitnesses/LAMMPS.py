import logging
import subprocess

import structopt
from structopt.tools import root, single_core, parallel
from structopt.tools.parallel import allgather


@parallel
def fitness(population):
    """Perform the LAMMPS fitness calculation on an entire population.

    Args:
        population (Population): the population to evaluate
    """
    if structopt.parameters.globals.USE_MPI4PY:
        logger = logging.getLogger('by-rank')
    else:
        logger = logging.getLogger('output')

    to_fit = [individual for individual in population if individual._modified]
    ncores = structopt.parameters.globals.ncores
    rank = structopt.parameters.globals.rank

    individuals_per_core = {r: [] for r in range(ncores)}
    for i, individual in enumerate(to_fit):
        individuals_per_core[i % ncores].append(individual.index)

    for index in individuals_per_core[rank]:
        individual = population[index]
        assert individual.index == index
        energy = individual.fitnesses.LAMMPS.get_energy(individual)
        individual.LAMMPS = energy
        logger.info('Individual {0} for LAMPPS evaluation had energy {1}'.format(i, energy))

    fits = [individual.LAMMPS for individual in population]

    if structopt.parameters.globals.ncores > 1:
        fits = allgather(fits, individuals_per_core)

    return fits

