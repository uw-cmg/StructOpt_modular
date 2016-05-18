import importlib

import structopt
from .crossovers import Crossovers
from .predators import Predators
from .selections import Selections
from .fitnesses import Fitnesses
from .relaxations import Relaxations
from .mutations import Mutations
from structopt.tools import root, single_core, parallel

from mpi4py import MPI


class Population(list):
    """A list-like class that contains the Individuals and the operations to be run on them."""

    @single_core
    def __init__(self):
        self.structure_type = structopt.parameters.generators.structure_type.lower()
        importlib.import_module('structopt.{}'.format(self.structure_type))
        self.crossovers = Crossovers()
        self.predators = Predators()
        self.selections = Selections()
        self.fitnesses = Fitnesses()
        self.relaxations = Relaxations()
        self.mutations = Mutations()

        # Import the structure type class: e.g from structopt.crystal import Crystal
        # Unfortunately `from` doesn't seem to work implicitly
        # so a getattr on the module is needed
        Structure = getattr(importlib.import_module('structopt.{}'.format(self.structure_type)),
                            self.structure_type.title())

        # Generate/load initial structures
        for structure_information in structopt.parameters.generators.initializers:
            for i in range(structure_information.number_of_individuals):
                structure = Structure(index=i, **structure_information.data)
                self.append(structure)

        self.total_number_of_individuals = len(self)


    def __getstate__(self):
        # Copy the object's state from self.__dict__ which contains
        # all our instance attributes. Always use the dict.copy()
        # method to avoid modifying the original state.
        state = self.__dict__.copy()
        # Remove the unpicklable entries.
        del state['structure_type']
        del state['crossovers']
        del state['predators']
        del state['selections']
        del state['fitnesses']
        del state['relaxations']
        del state['mutations']
        return state


    @single_core
    def replace(self, a_list):
        self.clear()
        self.extend(a_list)
        for i, individual in enumerate(self):
            individual.index = i


    @single_core
    def extend(self, other):
        super().extend(other)
        for i, individual in enumerate(self):
            individual.index = i

    @parallel
    def crossover(self):
        """Perform crossovers on the population."""
        children = self.crossovers.crossover(self)
        self.extend(children)
        self.crossovers.post_processing()


    @root
    def mutate(self):
        """Relax the entire population."""
        self.mutations.mutate(self)
        self.mutations.post_processing()
        return self


    @parallel
    def fitness(self):
        """Perform the fitness evaluations on the entire population."""
        fits = self.fitnesses.fitness(self)

        # Store the individuals total fitness for each individual
        for i, individual in enumerate(self):
            individual._fitness = fits[i]

        # Set each individual to unmodified so that the fitnesses wont't be recalculated
        for individual in self:
            individual._modifed = False

        self.fitnesses.post_processing()

        fits = MPI.COMM_WORLD.bcast(fits, root=0)
        return fits


    @parallel
    def relax(self):
        """Relax the entire population. """
        self.relaxations.relax(self)
        self.relaxations.post_processing()


    @root
    def kill(self, fits):
        """Remove individuals from the population whose fingerprints are very similar.
        The goal of this selection-like scheme is to encourage diversity in the population.

        Args:
            fits (list<float>): the fitnesses of each individual in the population
        """
        self.predators.select_predator()
        self.predators.kill(self, fits, nkeep=self.total_number_of_individuals)
        self.predators.post_processing()
        return self


    @root
    def select(self, fits):
        """Select the individuals in the population to keep for the next generation.

        Args:
            fits (list<float>): the fitnesses of each individual in the population
        """
        self.selections.select_selection()
        self.selections.select(self, fits)
        self.selections.post_processing()
        return self

