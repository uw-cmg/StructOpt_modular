import logging
import random

import structopt


class StructureType(object):
    """An abstract base class for a structure type."""

    def __init__(self, cls_name):
        self.logger = logging.getLogger('default')
        self.parameters = structopt.globals

        globals()['crossovers'] = import_module('structopt.{cls_name}.crossovers'.format(cls_name=cls_name))
        globals()['fingerprinters'] = import_module('structopt.{cls_name}.fingerprinters'.format(cls_name=cls_name))
        globals()['fitnesses'] = import_module('structopt.{cls_name}.fitnesses'.format(cls_name=cls_name))
        globals()['generators'] = import_module('structopt.{cls_name}.generators'.format(cls_name=cls_name))
        globals()['mutations'] = import_module('structopt.{cls_name}.mutations'.format(cls_name=cls_name))
        globals()['postprocessing'] = import_module('structopt.{cls_name}.postprocessing'.format(cls_name=cls_name))
        globals()['predators'] = import_module('structopt.{cls_name}.predators'.format(cls_name=cls_name))
        globals()['relaxations'] = import_module('structopt.{cls_name}.relaxations'.format(cls_name=cls_name))
        globals()['selections'] = import_module('structopt.{cls_name}.selections'.format(cls_name=cls_name))
        globals()['tools'] = import_module('structopt.{cls_name}.tools'.format(cls_name=cls_name))

        self.crossovers = crossovers.Crossovers(structopt.parameters.crossovers)
        self.fingerprinters = fingerprinters.Fingerprinters(structopt.parameters.fingerprinters)
        self.fitnesses = fitnesses.Fitnesses(structopt.parameters.fitnesses)
        self.generators = generators.Generators(structopt.parameters.generators)
        self.mutations = mutations.Mutations(structopt.parameters.mutations)
        self.postprocessing = postprocessing.Postprocessing(structopt.parameters.postprocessing)
        self.predators = predators.Predators(structopt.parameters.predators)
        self.relaxations = relaxations.Relaxations(structopt.parameters.relaxations)
        self.selections = selections.Selections(structopt.parameters.selections)
        self.tools = tools.Tools(structopt.parameters.tools)

        # Setup the output files

        # Set starting convergence
        self.converged = False

        # Prep output monitoring

        # Initialize random number seed
        random.seed(self.seed)

        # Write the input parameters to the output file
        self.logger.debug('Writing the input parameters to output file')
        structopt.fileio.parameters.write(self)

        # Allow structopt to see me
        structopt.parameters.globals.structuretype = self

    def crossover(self):
        self.crossovers.select_crossover()
        return self.crossovers.crossover(self.individuals)

    def mutate(self):
        self.mutations.select_mutation()
        return self.mutations.mutate(self.individuals)

    def relax(self):
        return self.relaxations.relax(self.individuals)

    def fitness(self):
        return self.fitness(self.individuals)

    def select(self):
        self.selections.select_selection()
        return self.selections.select(self.individuals)

    def kill(self):
        self.predators.select_predator()
        return self.predators.kill(self.individuals)

