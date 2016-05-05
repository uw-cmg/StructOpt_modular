import structopt.tools.structopt_lammps


class LAMMPS(object):
    def __init__(self):
        self.parameters = structopt.parameters.fitnesses.LAMMPS


    def get_command(self, individual):
        raise NotImplementedError


    def get_energy(self, individual):
        return structopt.tools.structopt_lammps.run(self.parameters, individual, relax=False)

