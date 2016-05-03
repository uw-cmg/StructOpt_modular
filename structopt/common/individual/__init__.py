import logging
import random
import ase
from importlib import import_module

import structopt


class Individual(ase.Atoms):
    """An abstract base class for a structure."""

    def __init__(self, **kwargs):
        # Load in the appropriate functionality
        fingerprinters = import_module('structopt.{cls_name}.fingerprinters'.format(cls_name=self.__class__.__name__.lower()))
        fitnesses = import_module('structopt.{cls_name}.fitnesses'.format(cls_name=self.__class__.__name__.lower()))
        mutations = import_module('structopt.{cls_name}.mutations'.format(cls_name=self.__class__.__name__.lower()))
        relaxations = import_module('structopt.{cls_name}.relaxations'.format(cls_name=self.__class__.__name__.lower()))

        self.fingerprinters = fingerprinters.Fingerprinters()
        self.fitnesses = fitnesses.Fitnesses()
        self.mutations = mutations.Mutations()
        self.relaxations = relaxations.Relaxations()

        # Initialize the ase.Atoms structure
        generators = import_module('structopt.{cls_name}.generators'.format(cls_name=self.__class__.__name__.lower()))
        generators.generate(self, **kwargs)

    def mutate(self):
        self.mutations.select_mutation()
        return self.mutations.mutate(self)

    def relax(self):
        return self.relaxations.relax(self)

    def fitness(self):
        return self.fitnesses.fitness(self)

    def fingerprint(self):
        return self.fingerprinters.fingerprint(self)

