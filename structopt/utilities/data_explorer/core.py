import os
import sys
import json
import re
import importlib
import weakref
from collections import defaultdict, Counter
import numpy as np

from common import lazy, lazyproperty

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
        self.output_file = os.path.join(dir, 'output.log')
        self._icache = {}
        self._created_killed = None

    @lazyproperty
    def generations(self):
        return Generations(self.genealogy_file, self)

    @lazyproperty
    def parameters(self):
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
        parameters = DictionaryObject(json.loads(parameters))
        sys.modules['gparameters'] = parameters
        return parameters

    @lazyproperty
    def _fitnesses(self):
        fitnesses = {}
        for line in open(self.fitnesses_file):
            firsthalf = re.compile(".+ : INFO : Generation (\d+), Individual (\d+):(.*)")
            secondhalf = re.compile("([\w]+): (\d+.\d+)")
            generation, id, fitness_str = firsthalf.findall(line)[0]
            fitnesses[(int(generation), int(id))] = {fit: float(value) for fit, value in secondhalf.findall(fitness_str)}
        return fitnesses

    def _get_fitness(self, id, generation):
        fits = self._fitnesses[(generation, id)]
        weights = {name: module.weight for name, module in self.parameters.fitnesses.items()}
        fitness = sum(fits[module]*weights[module] for module in fits)
        return fitness

    def _get_module_fitness(self, id, generation, module):
        return self._fitnesses[(generation, id)][module]

    def __getitem__(self, index):
        return self.generations[index]

    def __len__(self):
        return len(self.generations)

    def get_historical_individual(self, id):
        try:
            history = self._icache[id]
        except KeyError:
            history = IndividualHistory(id, self)
            self._icache[id] = history
        return history

    def _get_historical_info(self, id):
        if self._created_killed is None:
            #print("Parsing Genealogy file...")
            created_killed = defaultdict(lambda: [np.inf, -np.inf])
            for generation, population in enumerate(self):
                #sys.stdout.write("\r%d%%" % generation)
                #sys.stdout.flush()
                for id_ in population.keys():
                    if generation < created_killed[id_][0]:
                        created_killed[id_][0] = generation
                    if generation > created_killed[id_][1]:
                        created_killed[id_][1] = generation
            self._created_killed = created_killed
            #print('')
        return self._created_killed[id]


class Generations(object):
    def __init__(self, genealogy_filename, dataexplorer):
        self._dataexplorer = weakref.ref(dataexplorer)
        generations = open(genealogy_filename).readlines()
        self._data = [line.strip().split(' : INFO : ')[-1] for line in generations]

    def __getitem__(self, index):
        if isinstance(self._data[index], Population):
            return self._data[index]
        else:
            self._data[index] = Population(index, self._data[index], self._dataexplorer())
            return self._data[index]

    def __len__(self):
        return len(self._data)


class Population(dict):
    def __init__(self, generation, individuals_from_genealogy_file, dataexplorer):
        self._generation = generation
        self._dataexplorer = weakref.ref(dataexplorer)
        self._data = individuals_from_genealogy_file
        self._loaded = False

    @property
    def generation(self):
        return self._generation

    def __hash__(self):
        return self.generation

    def _load(self):
        self._data = self._data.split()
        self._data = [parse_gene(ind) for ind in self._data]
        self._data = {id: (c, m) for id, c, m in self._data}
        super().__init__({id: None for id in self._data if id not in self})
        self._loaded = True

    def __getitem__(self, key):
        if key not in self or super().__getitem__(key) is None:
            if not self._loaded:
                self._load()
                if key not in self._data:
                    raise KeyError(key)
            history = self._dataexplorer().get_historical_individual(key)
            ct, mt = self._data[key]
            self[key] = Individual(key, self.generation, history, self._dataexplorer(), ct, mt)
        return super().__getitem__(key)

    def __len__(self):
        return len(self._data)

    def keys(self):
        if not self._loaded:
            self._load()
        return super().keys()

    def items(self):
        if not self._loaded:
            self._load()
        for key in self:
            self[key]
        return super().items()

    def values(self):
        if not self._loaded:
            self._load()
        for key in self:
            self[key]
        return super().values()

    def __iter__(self):
        if not self._loaded:
            self._load()
        return super().__iter__()


class IndividualHistory(dict):
    def __init__(self, id, dataexplorer):
        self._id = id
        self._dataexplorer = weakref.ref(dataexplorer)
        self._created_on = None
        self._killed_on = None
        self._loaded = False

    def __repr__(self):
        return "<IndividualHistory {}>".format(self.id)
    __str__ = __repr__

    @property
    def id(self):
        return self._id

    def __hash__(self):
        return self.id

    def _load(self):
        super().__init__({generation: None for generation in range(self.created_on, self.killed_on+1)})

    def __getitem__(self, generation):
        return self._dataexplorer()[generation][self.id]

    @property
    def created_on(self):
        if self._created_on is not None:
            return self._created_on
        c, k = self._dataexplorer()._get_historical_info(self.id)
        self._created_on = c
        self._killed_on = k
        return self._created_on

    @property
    def killed_on(self):
        if self._killed_on is not None:
            return self._killed_on
        c, k = self._dataexplorer()._get_historical_info(self.id)
        self._created_on = c
        self._killed_on = k
        return self._killed_on

    def __contains__(self, item):
        if not self._loaded:
            self._load()
        return super().__contains__(item)

    def __iter__(self):
        if not self._loaded:
            self._load()
        return super().__iter__()


class Individual(object):
    def __init__(self, id, generation, history, dataexplorer, crossover_tag=None, mutation_tag=None):
        self._id = id
        self._history = history
        self._dataexplorer = weakref.ref(dataexplorer)
        self._generation = generation
        self._loaded = False

        self.mutation_tag = mutation_tag
        self.crossover_tag = crossover_tag
        if self.mutation_tag is not None:
            self.mutation_tag = self.mutation_tag
        if self.crossover_tag is not None:
            parents = re.compile('\((\d+)\+(\d+)\)(\w+)')
            parent1, parent2, tag = parents.findall(self.crossover_tag)[0]
            self.parent1 = int(parent1)
            self.parent2 = int(parent2)
            self.crossover_tag = tag

    def __hash__(self):
        return id(self)
        #return (self.id, self.generation)

    @property
    def id(self):
        return self._id
    @property
    def generation(self):
        return self._generation

    @property
    def history(self):
        return self._history

    @lazyproperty
    def fitness(self):
        return self._dataexplorer()._get_fitness(self.id, self.generation)

    @lazyproperty
    def LAMMPS(self):
        return self._dataexplorer()._get_module_fitness(self.id, self.generation, "LAMMPS")

    @lazyproperty
    def STEM(self):
        return self._dataexplorer()._get_module_fitness(self.id, self.generation, "STEM")

    @lazyproperty
    def FEMSIM(self):
        return self._dataexplorer()._get_module_fitness(self.id, self.generation, "FEMSIM")

    def __repr__(self):
        return "<Individual {} @{}>".format(self.id, self.generation)
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
        parameters = self._dataexplorer().parameters
        if "generators" in parameters:
            parameters.pop("generators")
        self.structure_type = parameters.structure_type.lower()
        module = importlib.import_module('structopt.{}'.format(self.structure_type))
        Structure = getattr(module, self.structure_type.title())
        filename = os.path.join(self._dataexplorer().XYZs, 'generation{}'.format(self.generation), 'individual{}.xyz'.format(self.id))
        generator_parameters = {"read_xyz": {"filename": filename}}
        self._structure = Structure(id=self.id,
                              relaxation_parameters=parameters.relaxations,
                              fitness_parameters=parameters.fitnesses,
                              mutation_parameters=parameters.mutations,
                              pso_moves_parameters=parameters.pso_moves,
                              generator_parameters=generator_parameters)
        self._loaded = True
        return self._structure

    def __getattr__(self, key):
        if not self._loaded:
            raise AttributeError("'{}' has no attribute '{}'".format(self.__class__.__name__, key))
        return getattr(self._structure, key)

    def get_parents(self, dataexplorer):
        mparent = None
        c1parent = None
        c2parent = None
        if self.crossover_tag is not None:
            c1parent = dataexplorer[self.generation-1][self.parent1]
            c2parent = dataexplorer[self.generation-1][self.parent2]
        elif self.mutation_tag is not None:
                mparent = dataexplorer[self.generation-1][self.id]
        return mparent, c1parent, c2parent

    def full_history(self, generations):
        if not hasattr(self, '_full_history'):
            history = []
            history.append((self.created_on, self))
            for parent in self.get_parents(generations):
                if parent is None: continue
                history.extend(parent.history(generations))
            self._full_history = history
        return self._full_history
        #history = defaultdict(list)
        #history[self.created_on].append(self)
        #for parent in self.get_parents(generations):
        #    if parent is None: continue
        #    history[parent.created_on].append(parent.history(generations))
        #return dict(history)

    def print_history(self, generations):
        history = self.full_history(generations)
        dhist = defaultdict(list)
        for gen, ind in history:
            dhist[gen].append(ind)
        for gen, hist in sorted(dhist.items()):
            print('{gen}: {hist}'.format(gen=gen, hist=', '.join(str(i) for i in set(hist))))

    def count_modifiers(self, generations):
        return self._count_modifiers(self.full_history(generations))

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

