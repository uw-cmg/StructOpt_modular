import logging
import subprocess

import structopt


def fitness(population):
    logger = logging.getLogger('output')
    chisqs = []
    if structopt.globals.parameters.rank == 0:
        # Setup each individual and get the command for each individual
        commands = []
        for individual in population:
            command = individual.fitnesses.femsim.get_command(individual)
            commands.append(command)

        # Run the parallelized `mpiexec` command
        commands = ['-np 1 {command}'.format(command=command) for command in commands]  # TODO: correctly allocate cores
        command = 'mpiexec {}'.format(' : '.join(commands))
        subprocess.call(command)

        # Collect the results for each chisq and return them
        for individual in population:
            chisq = individual.fitnesses.femsim.get_chisq(individual)
            logger.info('Individual {0} for FEMSIM evaluation had chisq {1}'.format(i, chisq))
            chisqs.append(chisq)
    return chisqs

