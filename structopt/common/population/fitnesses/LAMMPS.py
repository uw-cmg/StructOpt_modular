import logging
import subprocess

import structopt


def fitness(population):
    if structopt.parameters.globals.USE_MPI4PY:
        from mpi4py import MPI
        logger = logging.getLogger('by-rank')
    else:
        logger = logging.getLogger('output')

    energies = []
    for i, individual in enumerate(population):
        if structopt.parameters.globals.USE_MPI4PY and structopt.parameters.globals.rank == i:
            command = individual.fitnesses.LAMMPS.get_command(individual)
            subprocess.call(command, shell=True, stdout=subprocess.DEVNULL)
            energy = individual.fitnesses.LAMMPS.get_energy(individual)
            logger.info('Individual {0} for LAMPPS evaluation had energy {1}'.format(i, energy))
            energies = MPI.COMM_WORLD.gather(energy, root=0)
        else:
            command = individual.fitnesses.LAMMPS.get_command(individual)
            subprocess.call(command, shell=True, stdout=subprocess.DEVNULL)
            energy = individual.fitnesses.LAMMPS.get_energy(individual)
            logger.info('Individual {0} for LAMPPS evaluation had energy {1}'.format(i, energy))
            energies.append(energy)
 
    return energies

