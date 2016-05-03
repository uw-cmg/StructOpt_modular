import subprocess
from mpi4py import MPI

import structopt


def relax(population):
    for i, individual in enumerate(population):
        if structopt.parameters.globals.rank == i:
            print(individual.relaxations.LAMMPS)
            individual.relaxations.LAMMPS.relax(individual)
    return None

