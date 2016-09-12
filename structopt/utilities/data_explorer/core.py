import os
import json
import re
import importlib
from collections import defaultdict, Counter

from structopt.io import read_xyz
from structopt.tools.dictionaryobject import DictionaryObject


def parse_gene(gene):
    gene_tester = re.compile('(?P<cross>\(\d+\+\d+\)\D+)?(?P<id>\d+)(?P<mut>m\D+)?')
    crossover, id, mutation = gene_tester.findall(gene)[0]
    return int(id), crossover or None, mutation or None


class DataExplorer(object):
    def __init__(self, dir):
        self.genealogy_file = os.path.join(dir, 'genealogy.log')
        self.fitnesses_file = os.path.join(dir, 'fitnesses.log')
        self.XYZs = os.path.join(dir, 'XYZs')

        # Load input parameters from output file
        self.output_file = os.path.join(dir, 'output.log')
        parameters = []
        start_flag = False
        for line in open(self.output_file):
            if 'Current parameters:' in line:
                start_flag = True
                continue
            if start_flag:
                if 'INFO : {' in line:
                    line = line.split('INFO : ')[1]
                parameters.append(line)
            if line == '}\n':
                start_flag = False
                break
        parameters = ''.join(parameters)
        self.parameters = DictionaryObject(json.loads(parameters))
        import logging
        logging.parameters = self.parameters.logging

        self.generations = Generations(self.genealogy_file)
        for population in self.generations:
            for individual in population.values():
                individual.filename = os.path.join(self.XYZs, 'generation{}'.format(individual.created_on), 'individual{}.xyz'.format(individual.id))
                individual.parameters = DictionaryObject(self.parameters.copy())


class Generations(list):
    def __init__(self, genealogy_filename):
        generations = open(genealogy_filename).readlines()
        generations = [line.strip().split(' : INFO : ')[-1] for line in generations]
        super().__init__(line.split() for line in generations)
        individuals = defaultdict(bool)
        for i, population in enumerate(self):
            for j, gene in enumerate(population):
                id, crossover, mutation = parse_gene(gene)
                if mutation is None:
                    individuals[id] = individuals[id] or Individual(id=id, created_on=i, mutation_tag=mutation, crossover_tag=crossover)
                else:
                    individuals[id] = Individual(id=id, created_on=i, mutation_tag=mutation, crossover_tag=crossover)
                population[j] = individuals[id]
            self[i] = Population(self[i])
            for individual in list(individuals.keys()):
                if individual not in self[i]:
                    individuals[individual].killed_on = i
                    individuals.pop(individual)
        for individual in individuals:
            individuals[individual].killed_on = len(self)


class Population(dict):
    def __init__(self, individuals):
        super().__init__([(individual.id, individual) for individual in individuals])


class Individual(object):
    def __init__(self, id, created_on, killed_on=None, mutation_tag=None, crossover_tag=None, filename=None, parameters=None):
        self.id = id
        self.created_on = created_on
        self.killed_on = killed_on
        self.filename = filename
        self.mutation_tag = mutation_tag
        self.crossover_tag = crossover_tag
        if self.mutation_tag is not None:
            self.mutation_tag = self.mutation_tag[1:]
        if self.crossover_tag is not None:
            parents = re.compile('\((\d+)\+(\d+)\)(\w+)')
            parent1, parent2, tag = parents.findall(self.crossover_tag)[0]
            self.parent1 = int(parent1)
            self.parent2 = int(parent2)
            self.crossover_tag = tag
        self._loaded = False

    def __repr__(self):
        if self.mutation_tag is None:
            mtag = ''
        else:
            mtag = 'm' + self.mutation_tag
        if self.crossover_tag is None:
            return '{id}{mtag}'.format(id=self.id, mtag=mtag)
            #return '{id}{mtag} @{created_on}'.format(id=self.id, mtag=mtag, created_on=self.created_on)
        else:
            return '({p1}+{p2}){ctag}{id}{mtag}'.format(p1=self.parent1, p2=self.parent2, ctag=self.crossover_tag or '', id=self.id, mtag=mtag)
            #return '({p1}+{p2}){ctag}{id}{mtag} @{created_on}'.format(p1=self.parent1, p2=self.parent2, ctag=self.crossover_tag or '', id=self.id, mtag=mtag, created_on=self.created_on)
    __str__ = __repr__

    def load_structure(self):
        if "generators" in self.parameters:
            self.parameters.pop("generators")
        self.structure_type = self.parameters.structure_type.lower()
        module = importlib.import_module('structopt.{}'.format(self.structure_type))
        Structure = getattr(module, self.structure_type.title())
        generator_parameters = {"read_xyz": {"filename": self.filename}}
        self._structure = Structure(id=self.id,
                              relaxation_parameters=self.parameters.relaxations,
                              fitness_parameters=self.parameters.fitnesses,
                              mutation_parameters=self.parameters.mutations,
                              pso_moves_parameters=self.parameters.pso_moves,
                              generator_parameters=generator_parameters)
        self._loaded = True
        return self._structure

    def __getattr__(self, key):
        if not self._loaded:
            raise AttributeError("'{}' has no attribute '{}'".format(self.__class__.__name__, key))
        return getattr(self._structure, key)

    def get_parents(self, generations):
        mparent = None
        c1parent = None
        c2parent = None
        if self.crossover_tag is not None:
            c1parent = generations[self.created_on-1][self.parent1]
            c2parent = generations[self.created_on-1][self.parent2]
        elif self.mutation_tag is not None:
                mparent = generations[self.created_on-1][self.id]
        return mparent, c1parent, c2parent


    def history(self, generations):
        if not hasattr(self, '_history'):
            history = []
            history.append((self.created_on, self))
            for parent in self.get_parents(generations):
                if parent is None: continue
                history.extend(parent.history(generations))
            self._history = history
        return self._history
        #history = defaultdict(list)
        #history[self.created_on].append(self)
        #for parent in self.get_parents(generations):
        #    if parent is None: continue
        #    history[parent.created_on].append(parent.history(generations))
        #return dict(history)


    def print_history(self, generations):
        history = self.history(generations)
        dhist = defaultdict(list)
        for gen, ind in history:
            dhist[gen].append(ind)
        for gen, hist in sorted(dhist.items()):
            print('{gen}: {hist}'.format(gen=gen, hist=', '.join(str(i) for i in set(hist))))


    def count_modifiers(self, generations):
        return self._count_modifiers(self.history(generations))

    @staticmethod
    def _count_modifiers(history):
        crossovers = Counter()
        mutations = Counter()
        for gen, individual in history:
            crossovers[individual.crossover_tag] += 1
            mutations[individual.mutation_tag] += 1
        crossovers.pop(None, None)
        mutations.pop(None, None)
        return crossovers, mutations

