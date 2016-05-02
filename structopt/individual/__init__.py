import logging
import random
import import_module
import ase

import structopt.parameters


class Individual(ase.Atoms):
    """An abstract base class for a structure."""

    def __init__(self, *args, **kwargs):
        # Load in the appropriate functionality
        fingerprinters = import_module('structopt.{cls_name}.fingerprinters'.format(cls_name=self.__class__.__name__.lower()))
        fitnesses = import_module('structopt.{cls_name}.fitnesses'.format(cls_name=self.__class__.__name__.lower()))
        mutations = import_module('structopt.{cls_name}.mutations'.format(cls_name=self.__class__.__name__.lower()))
        relaxations = import_module('structopt.{cls_name}.relaxations'.format(cls_name=self.__class__.__name__.lower()))

        self.fingerprinters = fingerprinters.Fingerprinters(structopt.parameters.fingerprinters)
        self.fitnesses = fitnesses.Fitnesses(structopt.parameters.fitnesses)
        self.mutations = mutations.Mutations(structopt.parameters.mutations)
        self.relaxations = relaxations.Relaxations(structopt.parameters.relaxations)

        # Initialize the ase.Atoms structure
        generators = import_module('structopt.{cls_name}.generators'.format(cls_name=self.__class__.__name__))
        generators.generate(self, *args, **kwargs)

    def mutate(self):
        self.mutations.select_mutation()
        return self.mutations.mutate(self)

    def relax(self):
        return self.relaxations.relax(self)

    def fitness(self):
        return self.fitness(self)

    def fingerprint(self):
        return self.fingerprinters.fingerprint(self)

