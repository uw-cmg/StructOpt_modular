from structopt.tools import root, single_core, parallel
from structopt.tools.lammps import LAMMPS as lammps

class LAMMPS(object):
    """ """

    @single_core
    def __init__(self, parameters):
        # These variables never change
        self.parameters = parameters


    @single_core
    def get_command(self, individual):
        raise NotImplementedError


    @single_core
    def fitness(self, individual, generation=None):
        # Don't rerun lammps if:
        # 1) the individual is unmodified
        # 2) the energy has already been calculated via the relaxation
        if individual._relaxed and hasattr(individual, 'LAMMPS') and individual.LAMMPS is not None:
            return individual.LAMMPS
        else:
            print("Individual {} did not have an value for .LAMMPS or it was modified".format(individual.index))
            if generation is not None:
                calcdir = os.path.join(os.getcwd(), 'fitness-files/LAMMPS/generation-{}/individual-{}')
                calcdir = calcdir.format(generation, individual.index)
            else:
                calcdir = None

            calc = lammps(self.parameters, calcdir=calcdir)
            individual.set_calculator(calc)
            try:
                E = individual.get_potential_energy()
                print("Finished calculating fitness of individual {} on rank {} with LAMMPS".format(individual.index, rank))
            except RuntimeError:
                E = 0
                print("Error calculating fitness of individual {} on rank {} with LAMMPS".format(individual.index, rank))

            return E


