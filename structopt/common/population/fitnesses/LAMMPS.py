import logging
import subprocess
from mpi4py import MPI

import structopt


def fitness(population):
    logger = logging.getLogger('by-rank')
    for i, individual in enumerate(population):
        if structopt.parameters.globals.rank == i:
            command = individual.fitnesses.LAMMPS.get_command(individual)
            subprocess.call(command)
            energy = individual.fitnesses.LAMMPS.get_energy(individual)
            logger.info('Individual {0} for LAMPPS evaluation had energy {1}'.format(i, energy))
            energies = MPI.COMM_WORLD.allgather(energy)

    return energies

