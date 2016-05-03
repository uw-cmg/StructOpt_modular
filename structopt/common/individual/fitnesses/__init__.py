import structopt


class Fitnesses(object):
    def __init__(self):
        self.parameters = structopt.parameters.fitnesses
        self.modules = []

        if 'FEMSIM' in self.parameters.modules:
            from .FEMSIM import FEMSIM
            self.FEMSIM = FEMSIM()
            self.modules.append(self.FEMSIM)

        if 'LAMMPS' in self.parameters.modules:
            from .LAMMPS import LAMMPS
            self.LAMMPS = LAMMPS()
            self.modules.append(self.LAMMPS)


    def fitness(self, individual):
        fit = 0.0
        for i, module in enumerate(self.modules):  # TODO: ERROR: The implementation above will not sync self.modules with parameters.weights
            fit += module.fitness(individual) * self.parameters.weights[i]
        return fit

