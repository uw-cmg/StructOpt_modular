import structopt
from .LAMMPS import LAMMPS
#from .hard_sphere_cutoff import HardSphereCutoff
from .hard_sphere_cutoff import hard_sphere_cutoff
from structopt.tools import root, single_core, parallel


class Relaxations(object):
    """ """

    @single_core
    def __init__(self):
        self.parameters = structopt.parameters.relaxations
        self.modules = []

        for module in self.parameters.modules:
            setattr(self, module, globals()[module]())  # Initialized the class that was imported at the top of the file
            self.modules.append(getattr(self, module))


    @parallel
    def relax(self, individual):
        """Relax an individual.

        Args:
            individual (Individual): the individual to relax
        """
        for module in self.modules:
            module.relax(individual)
        individual._relaxed = True
        individual._fitted = False
        return None


    @single_core
    def post_processing(self):
        pass

