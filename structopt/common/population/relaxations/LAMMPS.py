import subprocess

import structopt


def relax(population):
    to_relax = [individual for individual in population if individual._modified]
    ncores = structopt.parameters.globals.ncores
    for i, individual in enumerate(to_relax):
        if structopt.parameters.globals.rank % structopt.parameters.globals.ncores == 0:
            individual.relaxations.LAMMPS.relax(individual)

