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
                structure = Structure(index=i, **structure_information.data)
                self.append(structure)

        self.total_number_of_individuals = len(self)


    def replace(self, a_list):
        self.clear()
        self.extend(a_list)


    def crossover(self):
        children = self.crossovers.crossover(self)
        self.extend(children)
        self.crossovers.post_processing()


    def select(self, fits):
        self.selections.select_selection()
        self.selections.select(self, fits, nkeep=self.total_number_of_individuals)
        self.selections.post_processing()


    def kill(self):
        self.predators.select_predator()
        self.predators.kill(self)
        self.predators.post_processing()


    def fitness(self):
        fits = self.fitnesses.fitness(self)
        for i, individual in enumerate(self):
            individual._fitness = fits[i]
        self.fitnesses.post_processing()
        return fits


    def relax(self):
        self.relaxations.relax(self)
        self.relaxations.post_processing()


    def mutate(self):
        self.mutations.mutate(self)
        self.mutations.post_processing()

