import subprocess

import structopt
from structopt.tools import root, single_core, parallel


@parallel
def relax(population):
    """Relax the entire population using a hard-sphere cutoff method.

    Args:
        population (Population): the population to relax
    """
    to_relax = [individual for individual in population if individual._modified]
    ncores = structopt.parameters.globals.ncores
    for i, individual in enumerate(to_relax):
        if structopt.parameters.globals.rank % structopt.parameters.globals.ncores == 0:
            individual.relaxations.LAMMPS.relax(individual)
    # TODO MPI collect

