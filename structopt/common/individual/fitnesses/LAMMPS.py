import numpy as np
from scipy.interpolate import interp1d
import os

from structopt.tools import single_core
from structopt.tools.lammps import LAMMPS as lammps
import gparameters


class LAMMPS(object):
    """ """

    @single_core
    def __init__(self, parameters):
        # These variables never change
        self.parameters = parameters
        self.output_dir = gparameters.logging.path

        # Set default normalization to E = E/natoms
        if "normalize" not in self.parameters:
            self.parameters.setdefault("normalize", {})
        self.parameters.normalize.setdefault('natoms', True)


    @single_core
    def get_command(self, individual):
        raise NotImplementedError


    @single_core
    def fitness(self, individual):
        # Don't rerun lammps if:
        # 1) the individual is unmodified
        # 2) the energy has already been calculated via the relaxation
        if individual._relaxed and hasattr(individual, 'LAMMPS') and individual.LAMMPS is not None:
            E = individual.LAMMPS
        else:
            print("Individual {} did not have an value for .LAMMPS or it was modified".format(individual.id))
            calcdir = os.path.join(self.output_dir, 'fitness/LAMMPS/generation{}/individual{}'.format(gparameters.generation, individual.id))
            rank = gparameters.mpi.rank

            calc = lammps(self.parameters, calcdir=calcdir)
            individual.set_calculator(calc)
            try:
                E = individual.get_potential_energy()
                print("Finished calculating fitness of individual {} on rank {} with LAMMPS".format(individual.id, rank))
            except RuntimeError:
                E = 0
                print("Error calculating fitness of individual {} on rank {} with LAMMPS".format(individual.id, rank))

        E = self.reference(E, individual)
        E = self.normalize(E, individual)
        individual.LAMMPS = E
        return E


    @single_core
    def reference(self, E, individual):
        """References the energy of the cluster to a reference energy"""
        if 'reference' not in self.parameters:
            return E

        E_refs = self.parameters.reference

        if type(self.parameters.reference) is float:
            return E - self.parameters.reference

        elif isinstance(self.parameters.reference, dict):
            symbols = individual.get_chemical_symbols()
            E_ref = [E_refs[symbol] for symbol in symbols]
            return E - sum(E_ref)

        else:
            raise TypeError('LAMMPS parameter "reference" must be either float or dict')

    @single_core
    def normalize(self, E, individual):
        if 'normalize' not in self.parameters:
            return E

        norms = self.parameters.normalize

        if 'natoms' in norms and norms['natoms'] is True:
            E /= len(individual)

        if 'particle' in norms:
            n = len(individual)
            f = np.poly1d(norms['particle'])
            E -= f(np.log(n))

        if 'alloy' in norms:
            element = norms['alloy']['element']
            xs = norms['alloy']['x']
            Es = norms['alloy']['E']

            linear = np.poly1d(np.polyfit([xs[0], xs[-1]], [Es[0], Es[-1]], 1))
            Es = np.array(Es) - linear(xs)

            f = interp1d(xs, Es)

            x = individual.get_chemical_symbols().count(element) / len(individual)
            E -= f(x)

        return E
