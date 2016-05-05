import logging
import subprocess

import structopt


def fitness(population):
    logger = logging.getLogger('output')
    chisqs = []
    if structopt.parameters.globals.rank == 0:
        # Setup each individual and get the command for each individual
        commands = []
        for individual in population:
            command = individual.fitnesses.FEMSIM.get_command(individual)
            commands.append(command)

        # Run the parallelized `mpiexec` command
        commands = ['-np 1 {command}'.format(command=command) for command in commands]  # TODO: correctly allocate cores
        command = 'mpiexec {}'.format(' : '.join(commands))
        subprocess.call(command, shell=True, stdout=subprocess.DEVNULL)

        # Collect the results for each chisq and return them
        for i, individual in enumerate(population):
            chisq = individual.fitnesses.FEMSIM.get_chisq(individual)
            logger.info('Individual {0} for FEMSIM evaluation had chisq {1}'.format(i, chisq))
            chisqs.append(chisq)
    return chisqs

