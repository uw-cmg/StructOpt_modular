import importlib
import numpy as np
from collections import Counter

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
    def __init__(self, parameters, individuals=None):
        self.parameters = parameters
        self.structure_type = self.parameters.generators.structure_type.lower()
        importlib.import_module('structopt.{}'.format(self.structure_type))
        self.crossovers = Crossovers(self.parameters.crossovers)
        self.predators = Predators(self.parameters.predators)
        self.selections = Selections(self.parameters.selections)
        self.fitnesses = Fitnesses(self.parameters.fitnesses)
        self.relaxations = Relaxations(self.parameters.relaxations)
        self.mutations = Mutations(self.parameters.mutations)
        self.generation = 0

        if individuals is None:
            # Import the structure type class: e.g from structopt.crystal import Crystal
            # Unfortunately `from` doesn't seem to work implicitly
            # so a getattr on the module is needed
            module = importlib.import_module('structopt.{}'.format(self.structure_type))
            Structure = getattr(module, self.structure_type.title())

            # Generate/load initial structures
            starting_index = 0
            for generator in self.parameters.generators.initializers:
                n = self.parameters.generators.initializers[generator].number_of_individuals
                for j in range(n):
                    kwargs = self.parameters.generators.initializers[generator].kwargs
                    generator_parameters = self.parameters.generators.initializers[generator]

                    # For the read_xyz, the input is a list of filenames. These need to be
                    # passed as arguments one by one instead of all at once
                    if generator == 'read_xyz':
                        kwargs = {'filename': kwargs[j]}

                    structure = Structure(index=starting_index + j,
                                          relaxation_parameters=self.parameters.relaxations,
                                          fitness_parameters=self.parameters.fitnesses,
                                          mutation_parameters=self.parameters.mutations,
                                          generator_parameters=generator_parameters)
                    self.append(structure)
                starting_index += n
        else:
            self.extend(individuals)

        self.initial_number_of_individuals = len(self)

    def __getstate__(self):
        # Copy the object's state from self.__dict__ which contains
        # all our instance attributes. Always use the dict.copy()
        # method to avoid modifying the original state.
        state = self.__dict__.copy()
        # Remove the unpicklable entries.
        del state['parameters']
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
        # TODO Make this call tool/parallel.allgather rather than reimplement it
        from mpi4py import MPI
        # The lists in individuals_per_core all need to be of the same length 
        max_individuals_per_core = max(len(individuals) for individuals in individuals_per_core.values())
        for rank, individuals in individuals_per_core.items():
            while len(individuals) < max_individuals_per_core:
                individuals.append(None)
 
        populations_per_rank = MPI.COMM_WORLD.allgather(self)
        #correct_population = [None for _ in range(sum(len(l) for l in populations_per_rank))]
        #correct_population = [None for _ in range(np.amax(list(individuals_per_core.values()))+1)]
        #correct_population = [None for _ in range(max_individuals_per_core+1)]
        correct_population = [None for _ in self]
        for rank, indices in individuals_per_core.items():
            for index in indices:
                if index is not None:
                    assert populations_per_rank[rank][index].index == index
                    correct_population[index] = populations_per_rank[rank][index]

        # If something didn't get sent, use the value on the core
        for i, individual in enumerate(correct_population):
            if individual is None:
                correct_population[i] = self[i]

        for i, individual in enumerate(correct_population):
            # See Individual.__getstate__; the below attributes don't get passed via MPI
            # so we need to reset them
            # TODO How can I make this dynamic? Can I call __getstate__???
            individual.fitnesses = self[i].fitnesses
            individual.relaxations = self[i].relaxations
            individual.mutations = self[i].mutations
            individual._calc = self[i]._calc
            individual._generator_kwargs = self[i]._generator_kwargs

        self.replace(correct_population)


    @single_core
    def replace(self, a_list):
        counter = Counter(individual.index for individual in self)
        assert max(counter.values()) == 1

        # TODO somehow check that this works
        mark_for_removal = [True for _ in self]
        for individual in a_list:
            for i, old_individual in enumerate(self):
                if individual.index == old_individual.index:
                    old_individual.update(individual)
                    mark_for_removal[i] = False
                    break

        mark_for_removal = [i for i, tf in enumerate(mark_for_removal) if tf]
        mark_for_removal.reverse()
        for i in mark_for_removal:
            self.pop(i)

        for i, individual in enumerate(self):
            individual.index = i


    @single_core
    def extend(self, other):
        super().extend(other)
        for i, individual in enumerate(self):
            individual.index = i

    @root
    def crossover(self, pairs):
        """Perform crossovers on the population."""
        children = self.crossovers.crossover(pairs)
        return children


    @root
    def mutate(self):
        """Perform mutations on the population."""
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
        """Relax the entire population."""
        self.relaxations.relax(self)


    @root
    def kill(self, fits):
        """Remove individuals from the population whose fingerprints are very similar.
        The goal of this selection-like scheme is to encourage diversity in the population.

        Args:
            fits (list<float>): the fitnesses of each individual in the population
        """
        self.predators.select_predator()
        self.predators.kill(self, fits, nkeep=self.initial_number_of_individuals)
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

