import structopt
import structopt.tools.structopt_lammps
from structopt.tools import root, single_core, parallel


class LAMMPS(object):
    """ """

    @single_core
    def __init__(self):
        # These variables never change
        self.parameters = structopt.parameters.relaxations.LAMMPS


    @single_core
    def get_command(self, individual):
        raise NotImplementedError


    @parallel
    def relax(self, individual):
        """Relax an individual.

        Args:
            individual (Individual): the individual to relax
        """
        print("Relaxing individual {} on rank {} with LAMMPS".format(individual.index, structopt.parameters.globals.rank))
        ret = structopt.tools.structopt_lammps.run(self.parameters, individual, relax=True)
        print("Finished relaxing individual {} on rank {} with LAMMPS".format(individual.index, structopt.parameters.globals.rank))
        return ret

