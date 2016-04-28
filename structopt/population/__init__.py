import import_module

import structopt.parameters
from structopt.population.crossovers import Crossovers
from structopt.population.predators import Predators
from structopt.population.selections import Selections


class Population(list):
    def __init__(self):
        self.crossovers = Crossovers(structopt.parameters.crossovers)
        self.predators = Predators(structopt.parameters.predators)
        self.selections = Selections(structopt.parameters.selections)

        # Generate/load initial structures
        for structure_information in structopt.parameters.generators:
            # Import the correct structure type class: e.g. structopt.crystal.Crystal
            Structure = import_module('structopt.{structure_type}.{title_case}'.format(structure_type=structure_information.structure_type.lower(),
                                                                                       title_case=structure_information.structure_type.title()))

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

