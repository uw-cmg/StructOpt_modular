import importlib

import structopt
from structopt.common.population.crossovers import Crossovers
from structopt.common.population.predators import Predators
from structopt.common.population.selections import Selections
from structopt.common.population.fitnesses import Fitnesses
from structopt.common.population.relaxations import Relaxations
from structopt.common.population.mutations import Mutations



class Population(list):
    """A list-like class that contains the Individuals and the operations to be run on them."""

    def __init__(self):
        self.structure_type = structopt.parameters.generators.structure_type.lower()
        importlib.import_module('structopt.{structure_type}'.format(structure_type=self.structure_type))
        self.crossovers = Crossovers()
        self.predators = Predators()
        self.selections = Selections()
        self.fitnesses = Fitnesses()
        self.relaxations = Relaxations()
        self.mutations = Mutations()

        # Import the correct structure type class: e.g. from structopt.crystal import Crystal
        # Unfortunately `from` doesn't seem to work implicitly so a getattr on the module is needed
        Structure = getattr(
            importlib.import_module('structopt.{structure_type}'.format(structure_type=self.structure_type)),
            self.structure_type.title()
        )
        # Generate/load initial structures
        for structure_information in structopt.parameters.generators.initializers:
            for i in range(structure_information.number_of_individuals):
                structure = Structure(**structure_information.data)
                self.append(structure)


    def crossover(self):
        self.crossovers.select_crossover()
        return self.crossovers.crossover(self)


    def select(self):
        self.selections.select_selection()
        return self.selections.select(self)


    def kill(self):
        self.predators.select_predator()
        return self.predators.kill(self)


    def fitness(self):
        return self.fitnesses.fitness(self)


    def relax(self):
        return self.relaxations.relax(self)


    def mutate(self):
        return self.mutations.mutate(self)

