import structopt
import structopt.tools.structopt_lammps
from structopt.tools import root, single_core, parallel


class LAMMPS(object):
    """ """

    @single_core
    def __init__(self, parameters=None):
        # These variables never change
        self.parameters = parameters or structopt.parameters.relaxations.LAMMPS


    @single_core
    def get_command(self, individual):
        raise NotImplementedError


    @parallel
    def relax(self, individual):
        """Relax an individual.

        Args:
            individual (Individual): the individual to relax
        """
        try:
            rank = structopt.parameters.globals.rank
        except:
            rank = 0
        print("Relaxing individual {} on rank {} with LAMMPS".format(individual.index, rank))
        ret = structopt.tools.structopt_lammps.run(self.parameters, individual, relax=True)
        print("Finished relaxing individual {} on rank {} with LAMMPS".format(individual.index, rank))
        return ret

