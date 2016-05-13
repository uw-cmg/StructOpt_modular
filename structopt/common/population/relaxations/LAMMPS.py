import subprocess

import structopt


@parallel
def relax(population):
    to_relax = [individual for individual in population if individual._modified]
    ncores = structopt.parameters.globals.ncores
    # TODO MPI send the individuals out to their respective cores
    for i, individual in enumerate(to_relax):
        if structopt.parameters.globals.rank % structopt.parameters.globals.ncores == 0:
            individual.relaxations.LAMMPS.relax(individual)
    # TODO MPI collect

