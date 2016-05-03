import importlib

import structopt
import structopt.cluster
from structopt.population.crossovers import Crossovers
from structopt.population.predators import Predators
from structopt.population.selections import Selections
import structopt.population.fitnesses



class Population(list):
    """A list-like class that contains the Individuals and the operations to be run on them."""

    def __init__(self):
        self.structure_type = structopt.parameters.generators.structure_type.lower()
        self.crossovers = Crossovers()
        self.predators = Predators()
        self.selections = Selections()
        self.fitnesses = [
            getattr(
                importlib.import_module('structopt.{structure_type}.fitnesses'.format(structure_type=self.structure_type)),
                module
            )
            for module in structopt.parameters.fitnesses.modules
        ]

        # Generate/load initial structures
        Structure = getattr(
            importlib.import_module('structopt.{structure_type}'.format(structure_type=self.structure_type)),
            self.structure_type.title()
        )
        for structure_information in structopt.parameters.generators.initializers:
            # Import the correct structure type class: e.g. from structopt.crystal import Crystal
            # Unfortunately `from` doesn't seem to work implicitly so a getattr on the module is needed
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
        fitnesses = np.zeros((len(self),), dtype=np.float)
        for i, module in enumerate(self.fitnesses):
            fits = module.fitness(self)

            # Save the fitness value for each module to each individual
            for j, individual in enumerate(self):
                setattr(individual, module, fits[j])

            fits = np.multiply(fits, structopt.parameters.fitnesses.weights)
            fitnesses = np.add(fitnesses, fits)

        return fitnesses


    def relax(self):
        for individual in self:
            individual.relax()


    def mutate(self):
        for individual in self:
            individual.mutate()

