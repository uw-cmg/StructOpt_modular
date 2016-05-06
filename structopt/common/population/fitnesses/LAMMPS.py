import logging
import subprocess

import structopt


def fitness(population):
    to_fit = [individual for individual in population if individual._modified]

    if structopt.parameters.globals.USE_MPI4PY:
        from mpi4py import MPI
        logger = logging.getLogger('by-rank')
    else:
        logger = logging.getLogger('output')

    for i, individual in enumerate(to_fit):
        if structopt.parameters.globals.USE_MPI4PY and structopt.parameters.globals.rank % structopt.parameters.globals.ncores == 0:
            energy = individual.fitnesses.LAMMPS.get_energy(individual)
            individual.LAMMPS = energy
            logger.info('Individual {0} for LAMPPS evaluation had energy {1}'.format(i, energy))
        else:
            energy = individual.fitnesses.LAMMPS.get_energy(individual)
            individual.LAMMPS = energy
            logger.info('Individual {0} for LAMPPS evaluation had energy {1}'.format(i, energy))
 
    return [individual.LAMMPS for individual in population]

