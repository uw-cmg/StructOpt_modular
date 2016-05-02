import importlib

import structopt
import structopt.cluster
from structopt.population.crossovers import Crossovers
from structopt.population.predators import Predators
from structopt.population.selections import Selections
import structopt.population.fitnesses



class Population(list):
    def __init__(self):
        self.crossovers = Crossovers()
        self.predators = Predators()
        self.selections = Selections()

        # Generate/load initial structures
        print(structopt.parameters.generators[0])
        print(type(structopt.parameters.generators[0]))
        for structure_information in structopt.parameters.generators:
            # Import the correct structure type class: e.g. from structopt.crystal import Crystal
            # Unfortunately `from` doesn't seem to work implicitly so a getattr on the module is needed
            print(structure_information)
            Structure = getattr(
                importlib.import_module('structopt.{structure_type}'.format(structure_type=structure_information.structure_type.lower())),
                structure_information.structure_type.title()
            )

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
        for individual in self:
            individual.fitness()

    def relax(self):
        for individual in self:
            individual.relax()

    def mutate(self):
        for individual in self:
            individual.mutate()

