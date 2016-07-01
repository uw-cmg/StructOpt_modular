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

    to_fit = [individual for individual in population if not individual._fitted]
    if structopt.parameters.globals.USE_MPI4PY:
        ncores = structopt.parameters.globals.ncores
    else:
        ncores = 1
    rank = structopt.parameters.globals.rank

    individuals_per_core = {r: [] for r in range(ncores)}
    for i, individual in enumerate(to_fit):
        individuals_per_core[i % ncores].append(individual.index)

    for index in individuals_per_core[rank]:
        individual = population[index]
        assert individual.index == index
        print("Running LAMMPS fitness evaluation on individual {}".format(individual.index))
        energy = individual.fitnesses.LAMMPS.get_energy(individual)
        individual.LAMMPS = energy
        logger.info('Individual {0} after LAMPPS evaluation has energy {1}'.format(individual.index, energy))

    fits = [individual.LAMMPS for individual in population]
    if structopt.parameters.globals.USE_MPI4PY:
        fits = allgather(fits, individuals_per_core)

    # Save the fitness value for the module to each individual after they have been allgathered
    for i, fit in enumerate(fits):
        to_fit[i].LAMMPS = fit

    return [individual.LAMMPS for individual in population]

