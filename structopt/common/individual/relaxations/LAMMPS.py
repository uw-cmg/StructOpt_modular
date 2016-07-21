import os
import logging
import structopt

from structopt.tools.lammps import LAMMPS as lammps
from structopt.tools import root, single_core, parallel

class LAMMPS(object):
    """ """

    @single_core
    def __init__(self, parameters):
        # These variables never change
        self.parameters = parameters


    @single_core
    def get_command(self, individual):
        raise NotImplementedError


    @parallel
    def relax(self, individual, **kwargs):
        """Relax an individual.

        Args:
            individual (Individual): the individual to relax
        """
        rank = logging.parameters.rank
        print("Relaxing individual {} on rank {} with LAMMPS".format(individual.index, rank))

        if 'generation' in kwargs:
            calcdir = os.path.join(os.getcwd(), 'relaxation-files/LAMMPS/generation-{}/individual-{}')
            calcdir = calcdir.format(kwargs['generation'], individual.index)
        else:
            calcdir = None

        calc = lammps(self.parameters, calcdir=calcdir)
        individual.set_calculator(calc)
        try:
            E = individual.get_potential_energy()
            print("Finished relaxing individual {} on rank {} with LAMMPS".format(individual.index, rank))
        except RuntimeError:
            E = np.nan
            print("Error relaxing individual {} on rank {} with LAMMPS".format(individual.index, rank))

        individual.set_calculator()
        individual.LAMMPS = E

        return

