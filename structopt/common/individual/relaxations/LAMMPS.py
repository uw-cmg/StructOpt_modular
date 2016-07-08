import logging
import structopt
import structopt.tools.structopt_lammps
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
    def relax(self, individual):
        """Relax an individual.

        Args:
            individual (Individual): the individual to relax
        """
        rank = logging.parameters.rank
        print("Relaxing individual {} on rank {} with LAMMPS".format(individual.index, rank))
        ret = structopt.tools.structopt_lammps.run(self.parameters, individual, relax=True, use_mpi4py=self.parameters.use_mpi4py)
        print("Finished relaxing individual {} on rank {} with LAMMPS".format(individual.index, rank))
        return ret

