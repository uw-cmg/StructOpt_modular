import structopt.tools.structopt_lammps
from structopt.tools import root, single_core, parallel


class LAMMPS(object):
    """ """

    @single_core
    def __init__(self):
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
        return structopt.tools.structopt_lammps.run(self.parameters, individual, relax=True)

