import structopt.tools.structopt_lammps
from structopt.tools import root, single_core, parallel


class LAMMPS(object):
    """ """

    @single_core
    def __init__(self):
        self.parameters = structopt.parameters.fitnesses.LAMMPS


    @single_core
    def get_command(self, individual):
        raise NotImplementedError


    @single_core
    def get_energy(self, individual):
        # Don't rerun lammps if:
        # 1) the individual is unmodified
        # 2) the energy has already been calculated via the relaxation
        if individual._modified or 'LAMMPS' not in structopt.parameters.relaxations.modules:
            return structopt.tools.structopt_lammps.run(self.parameters, individual, relax=False)
        else:
            return individual.LAMMPS

