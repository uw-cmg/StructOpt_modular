import os
import logging

from structopt.tools import root, single_core, parallel
from structopt.tools.lammps import LAMMPS as lammps
from structopt.tools.dictionaryobject import DictionaryObject

class LAMMPS(object):
    """ """

    @single_core
    def __init__(self, parameters):
        # These variables never change
        self.parameters = parameters
        self.output_dir = logging.parameters.path

        # Set default normalization to E = E/natoms
        if "normalize" not in self.parameters:
            self.parameters.setdefault("normalize", DictionaryObject({}))
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
            if hasattr(logging, 'parameters'):
                calcdir = os.path.join(self.output_dir, 'fitness/LAMMPS/generation{}/individual{}'.format(logging.parameters.generation, individual.id))
            else:
                calcdir = None

            calc = lammps(self.parameters, calcdir=calcdir)
            individual.set_calculator(calc)
            try:
                E = individual.get_potential_energy()
                print("Finished calculating fitness of individual {} on rank {} with LAMMPS".format(individual.id, logging.parameters.rank))
            except RuntimeError:
                E = 0
                print("Error calculating fitness of individual {} on rank {} with LAMMPS".format(individual.id, logging.parameters.rank))

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

        elif type(self.parameters.reference) in [dict, DictionaryObject]:
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

        return E
