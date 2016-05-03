import structopt


class Relaxations(object):
    def __init__(self):
        self.parameters = structopt.parameters.relaxations
        self.modules = []

        if 'LAMMPS' in self.parameters.modules:
            from .LAMMPS import LAMMPS
            self.LAMMPS = LAMMPS()
            self.modules.append(self.LAMMPS )


    def relax(self, individual):
        for module in self.modules:
            module.relax(individual)
        return None

