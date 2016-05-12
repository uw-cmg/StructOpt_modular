import structopt
from .FEMSIM import FEMSIM
from .LAMMPS import LAMMPS


class Fitnesses(object):
    def __init__(self):
        self.parameters = structopt.parameters.fitnesses
        self.modules = []

        for module in self.parameters.modules:
            setattr(self, module, globals()[module]())  # Initialized the class that was imported at the top of the file
            self.modules.append(getattr(self, module))


    def fitness(self, individual):
        fit = 0.0
        for i, module in enumerate(self.modules):
            fit += module.fitness(individual) * self.parameters.weights[i]
        return fit

    def post_processing(self):
        pass

