import structopt
from .FEMSIM import FEMSIM
from .LAMMPS import LAMMPS
from structopt.tools import root, single_core, parallel


class Fitnesses(object):
    """ """

    @single_core
    def __init__(self, parameters):
        self.parameters = parameters
        self.modules = []

        for module in self.parameters.modules:
            # Initialize the class that was imported at the top of the file and append it to the modules list
            parameters = getattr(self.parameters, module)
            setattr(self, module, globals()[module](parameters=parameters))
            self.modules.append(getattr(self, module))


    @parallel
    def fitness(self, individual, **kwargs):
        """Perform the fitness calculations on an individual.

        Args:
            individual (Individual): the individual to evaluate
        """
        fit = 0.0
        for i, module in enumerate(self.modules):
            fit += module.fitness(individual, **kwargs) * self.parameters.weights[i]
        return fit


    @single_core
    def post_processing(self):
        pass

