import importlib

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


    @single_core
    def replace(self, a_list):
        self.clear()
        self.extend(a_list)


    @parallel
    def crossover(self):
        children = self.crossovers.crossover(self)
        self.extend(children)
        self.crossovers.post_processing()


    @root
    def mutate(self):
        self.mutations.mutate(self)
        self.mutations.post_processing()


    @parallel
    def fitness(self):
        fits = self.fitnesses.fitness(self)
        for i, individual in enumerate(self):
            individual._fitness = fits[i]
        self.fitnesses.post_processing()
        return fits


    @parallel
    def relax(self):
        self.relaxations.relax(self)
        self.relaxations.post_processing()


    @root
    def kill(self, fits):
        self.predators.select_predator()
        self.predators.kill(self, fits)
        self.predators.post_processing()


    @root
    def select(self, fits):
        self.selections.select_selection()
        self.selections.select(self, fits, nkeep=self.total_number_of_individuals)
        self.selections.post_processing()

