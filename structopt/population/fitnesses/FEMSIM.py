import subprocess

import structopt.parameters


def fitness(population):
    logger = logging.getLogger('output')
    chisqs = []
    if structopt.parameters.rank == 0:
        # Setup each individual and get the command for each individual
        commands = []
        for individual in population:
            command = individual.femsim.get_command(
                # TODO
            )
            commands.append(command)

        # Run the parallelized `mpiexec` command
        commands = ['-np 1 {command}'.format(command=command) for command in commands]  # TODO correctly allocate cores
        command = 'mpiexec {}'.format(' : '.join(commands))
        subprocess.call(command)

        # Collect the results for each chisq and return them
        for individual in population:
            chisq = individual.femsim.collect_chisq(
                # TODO
            )
            logger.info('Individual {0} for FEMSIM evaluation had chisq {1}'.format(i, chisq))
            chisqs.append(chisq)
    return chisqs
