import subprocess

import structopt


def relax(population):
    # TODO update this so that it correctly distributes which core relaxes which individual
    for i, individual in enumerate(population):
        if structopt.parameters.globals.rank == i:
            individual.relaxations.LAMMPS.relax(individual)

