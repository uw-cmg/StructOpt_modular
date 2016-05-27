import importlib
import numpy as np

import structopt
from .crossovers import Crossovers
from .predators import Predators
from .selections import Selections
from .fitnesses import Fitnesses
from .relaxations import Relaxations
from .mutations import Mutations
from structopt.tools import root, single_core, parallel


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


    @parallel
    def allgather(self, individuals_per_core):
        """Performs an MPI.allgather on self (the population) and updates the
        correct individuals that have been modified based on the inputs from 
        individuals_per_core.

        See stuctopt.tools.parallel.allgather for a similar function.
        """
        if structopt.parameters.globals.USE_MPI4PY:
            from mpi4py import MPI
            populations_per_rank = MPI.COMM_WORLD.allgather(self)
            #correct_population = [None for _ in range(sum(len(l) for l in populations_per_rank))]
            correct_population = [None for _ in range(np.amax(list(individuals_per_core.values()))+1)]
            for rank, indices in individuals_per_core.items():
                for index in indices:
                    assert populations_per_rank[rank][index].index == index
                    correct_population[index] = populations_per_rank[rank][index]

            # If something didn't get sent, use the value on the core
            for i, individual in enumerate(correct_population):
                if individual is None:
                    correct_population[i] = self[i]

            for i, individual in enumerate(correct_population):
                # See Individual.__getstate__; the below attributes don't get passed via MPI
                # so we need to reset them
                individual.fitnesses = self[i].fitnesses
                individual.relaxations = self[i].relaxations
                individual.mutations = self[i].mutations
                individual._calc = self[i]._calc
                individual._kwargs = self[i]._kwargs

            self.replace(correct_population)


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
    def crossover(self, pairs):
        """Perform crossovers on the population."""
        children = self.crossovers.crossover(pairs)
        self.extend(children)


    @root
    def mutate(self):
        """Relax the entire population."""
        self.mutations.mutate(self)
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
            individual._fitted = True

        return fits


    @parallel
    def relax(self):
        """Relax the entire population. """
        self.relaxations.relax(self)


    @root
    def kill(self, fits):
        """Remove individuals from the population whose fingerprints are very similar.
        The goal of this selection-like scheme is to encourage diversity in the population.

        Args:
            fits (list<float>): the fitnesses of each individual in the population
        """
        self.predators.select_predator()
        self.predators.kill(self, fits, nkeep=self.total_number_of_individuals)
        return self


    @root
    def select(self, fits):
        """Select the individuals in the population to perform crossovers on.

        Args:
            fits (list<float>): the fitnesses of each individual in the population
        """
        self.selections.select_selection()
        pairs = self.selections.select(self, fits)
        return pairs

