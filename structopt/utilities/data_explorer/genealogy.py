import re
from collections import Counter


class Individual(object):
    def __init__(self, id=None, mutation_tag=None, crossover_tag=None, filename=None):
        self.id = id
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

    def __repr__(self):
        if self.mutation_tag is None:
            mtag = ''
        else:
            mtag = 'm' + self.mutation_tag
        if self.crossover_tag is None:
            return '{id}{mtag}'.format(id=self.id, mtag=mtag)
        else:
            return '({p1}+{p2}){ctag}{id}{mtag}'.format(p1=self.parent1, p2=self.parent2, ctag=self.crossover_tag or '', id=self.id, mtag=mtag)
    __str__ = __repr__


def parse_gene(gene):
    gene_tester = re.compile('(?P<cross>\(\d+\+\d+\)\D+)?(?P<id>\d+)(?P<mut>m\D+)?')
    crossover, id, mutation = gene_tester.findall(gene)[0]
    return int(id), crossover or None, mutation or None


def load_genealogy(filename):
    generations = open(filename).readlines()
    generations = [line.strip().split(' : INFO : ')[-1] for line in generations]
    generations = [line.split() for line in generations]
    for population in generations:
        for i, gene in enumerate(population):
            id, crossover, mutation = parse_gene(gene)
            individual = Individual(id=id, mutation_tag=mutation, crossover_tag=crossover)
            population[i] = individual
    return generations


def trace_history(generations, id):
    history = [None for _ in generations]
    # Traverse the generations in reverse order
    for generation in range(len(generations)-1, 0, -1):
        population = generations[generation]
        for individual in population:
            if individual.id == id:
                history[generation] = individual
                if individual.crossover_tag is not None:
                    history1 = trace_history(generations, individual.parent1)
                    history2 = trace_history(generations, individual.parent2)
                    history1 = history1[:generation]
                    history2 = history2[:generation]
                    parent_history = zip(history1, history2)
                    for generation_, history_ in enumerate(parent_history):
                        history[generation_] = history_
                break
    return history


def _count_modifiers(maybe_iterable):
    mutations = Counter()
    crossovers = Counter()
    if not isinstance(maybe_iterable, Individual):# and hasattr(maybe_iterable, '__iter__'):
        for thing in maybe_iterable:
            if thing is None: continue
            if not isinstance(thing, Individual):
                crossovers_, mutations_ = _count_modifiers(thing)
                #print("Is iterable:", maybe_iterable, crossovers_, mutations_)
                for x in crossovers_:
                    crossovers[x] += crossovers_[x]
                for m in mutations_:
                    mutations[m] += mutations_[m]
            else:
                mutations[thing.mutation_tag] += 1
                crossovers[thing.crossover_tag] += 1
    else:
        thing = maybe_iterable
        if thing is not None:
            mutations[thing.mutation_tag] += 1
            crossovers[thing.crossover_tag] += 1
        #print("Not iterable:", maybe_iterable, crossovers, mutations)
    #print("Returning:", crossovers, mutations)
    return crossovers, mutations

def count_modifiers(history):
    mutations = Counter()
    crossovers = Counter()
    for generation, things in enumerate(history):
        crossovers_, mutations_ = _count_modifiers(things)
        for x in crossovers_:
            crossovers[x] += crossovers_[x]
        for m in mutations_:
            mutations[m] += mutations_[m]
    return crossovers, mutations


if __name__ == '__main__':
    generations = load_genealogy('genealogy.log')
    history = trace_history(generations, 10)
    for generation, history_ in enumerate(history):
        print('{}  {}'.format(generation, history_))
    crossovers, mutations = count_modifiers(history)
    print("Crossovers:", crossovers)
    print("Mutations:", mutations)

