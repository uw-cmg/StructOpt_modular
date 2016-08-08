import logging

from structopt.tools import root, single_core, parallel

@parallel
def fitness(population, parameters):
    """Performs the STEM fitness calculation on an entire population.

    Parameters
    ----------
        population : structopt.Population
    """

    if parameters.use_mpi4py:
        logger = logging.getLogger('by-rank')
        ncores = logging.parameters.ncores
    else:
        logger = logging.getLogger('output')
        ncores = 1

    rank = logging.parameters.rank

    individuals_per_core = {r: [] for r in range(ncores)}
    for i, individual in enumerate(population):
        individuals_per_core[i % ncores].append(individual)

    for individual in individuals_per_core[rank]:
        print("Running STEM fitness evaluation on individual {}".format(individual.id))
        chi2 = individual.fitnesses.STEM.fitness(individual)
        individual.STEM = chi2
        logger.info('Individual {0} after STEM evaluation has chi^2 {1}'.format(individual.id, chi2))

    fits = [individual.STEM for individual in population]
    positions_per_core = {rank: [population.position(individual) for individual in individuals] for rank, individuals in individuals_per_core.items()}

    if parameters.use_mpi4py:
        fits = allgather(fits, positions_per_core)

    # Save the fitness value for the module to each individual after they have been allgathered
    for i, fit in enumerate(fits):
        population.get_by_position(i).STEM = fit

    return [individual.STEM for individual in population]    
