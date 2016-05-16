import logging
import subprocess
import math

import structopt
from structopt.tools import root, single_core, parallel


@root
def fitness(population):
    """Perform the FEMSIM fitness calculation on an entire population.

    Args:
        population (Population): the population to evaluate
    """
    to_fit = [individual for individual in population if individual._modified]
    cores_per_individual = structopt.parameters.globals.ncores // len(to_fit)
    # Round cores_per_individual down to nearest power of 2
    if cores_per_individual == 0:
        # We have too more individuals than cores, so each fitness scheme needs to be run multiple times
        cores_per_individual = 1

    pow(2.0, math.floor(math.log2(cores_per_individual)))

    logger = logging.getLogger('output')
    if structopt.parameters.globals.rank == 0:
        # Setup each individual and get the command for each individual
        commands = []
        for individual in to_fit:
            command = individual.fitnesses.FEMSIM.get_command(individual)
            commands.append(command)

        # Run the parallelized `mpiexec` command
        commands = ['-np {cores} {command}'.format(command=command, cores=cores_per_individual) for command in commands]
        command = 'mpiexec {}'.format(' : '.join(commands))
        subprocess.call(command, shell=True, stdout=subprocess.DEVNULL)

        # Collect the results for each chisq and return them
        for i, individual in enumerate(population):
            print("RUNNING FEMSIM FOR INDIVIDUAL", individual)
            chisq = individual.fitnesses.FEMSIM.get_chisq(individual)
            individual.FEMSIM = chisq
            logger.info('Individual {0} for FEMSIM evaluation had chisq {1}'.format(i, chisq))
    return [individual.FEMSIM for individual in population]

