import logging
import random
import ase
from importlib import import_module
import numpy as np

import structopt


class Individual(ase.Atoms):
    """An abstract base class for a structure."""

    def __init__(self, index, **kwargs):
        self.index = index

        cls_name = self.__class__.__name__.lower()
        # Load in the appropriate functionality
        fingerprinters = import_module('structopt.{cls_name}.individual.fingerprinters'.format(cls_name=cls_name))
        fitnesses = import_module('structopt.{cls_name}.individual.fitnesses'.format(cls_name=cls_name))
        mutations = import_module('structopt.{cls_name}.individual.mutations'.format(cls_name=cls_name))
        relaxations = import_module('structopt.{cls_name}.individual.relaxations'.format(cls_name=cls_name))

        self.fingerprinters = fingerprinters.Fingerprinters()
        self.fitnesses = fitnesses.Fitnesses()
        self.mutations = mutations.Mutations()
        self.relaxations = relaxations.Relaxations()

        # Initialize the ase.Atoms structure
        super().__init__()
        generators = import_module('structopt.{cls_name}.individual.generators'.format(cls_name=cls_name))
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


    def get_atom_indices_within_distance_of_atom(self, atom_index, distance):
        dists = self.get_distances(atom_index, slice(None, None, None))
        return np.where(dists < distance)]

    def get_nearest_atom_indices(self, atom_index, count):
        dists = self.get_distances(atom_index, slice(None, None, None))
        return np.argsort(dists)[:count]

