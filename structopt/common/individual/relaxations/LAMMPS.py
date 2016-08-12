import os
import logging
import structopt
import numpy as np

from structopt.tools.lammps import LAMMPS as lammps
from structopt.tools import root, single_core, parallel
from structopt.common.individual.mutations.move_surface_atoms import move_surface_atoms

class LAMMPS(object):
    """ """

    @single_core
    def __init__(self, parameters):
        # These variables never change
        self.parameters = parameters
        if hasattr(logging, 'parameters'):
            self.output_dir = logging.parameters.path
        else:
            self.output_dir = os.getcwd()


    @single_core
    def get_command(self, individual):
        raise NotImplementedError


    @parallel
    def relax(self, individual):
        """Relax an individual.

        Args:
            individual (Individual): the individual to relax
        """
        rank = logging.parameters.rank
        print("Relaxing individual {} on rank {} with LAMMPS".format(individual.id, rank))

        if hasattr(logging, 'parameters'):
            calcdir = os.path.join(self.output_dir, 'relaxation/LAMMPS/generation{}/individual{}'.format(logging.parameters.generation, individual.id))
        else:
            calcdir = None

        calc = lammps(self.parameters, calcdir=calcdir)
        individual.set_calculator(calc)
        try:
            E = individual.get_potential_energy()
            print("Finished relaxing individual {} on rank {} with LAMMPS".format(individual.id, rank))
        except RuntimeError:
            E = np.nan
            print("Error relaxing individual {} on rank {} with LAMMPS".format(individual.id, rank))

        individual.LAMMPS = E

        if 'repair' in self.parameters and self.parameters['repair']:
            E = self.repair(individual, logging.parameters.generation)
            if E is not None:
                individual.LAMMPS = E

        return

    @parallel
    def repair(self, individual, generation):
        """Repairs an individual. Currently takes isolated atoms moves them next to
        a non-isolated atom"""

        rank = logging.parameters.rank

        n_moves = move_surface_atoms(individual, max_natoms=1.0, move_CN=3)
        if n_moves == 0:
            return

        if generation is not None:
            calcdir = os.path.join(os.getcwd(), '{}/relaxation/LAMMPS/generation{}/individual{}-repair')
            calcdir = calcdir.format(self.output_dir, generation, individual.id)
        else:
            calcdir = None

        calc = lammps(self.parameters, calcdir=calcdir)
        individual.set_calculator(calc)

        try:
            E = individual.get_potential_energy()
            print("Finished repairing individual {} on rank {} with LAMMPS".format(individual.id, rank))
        except RuntimeError:
            E = np.nan
            print("Error repairing individual {} on rank {} with LAMMPS".format(individual.id, rank))

        return E
