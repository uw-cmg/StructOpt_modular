import functools
import random
from itertools import accumulate, combinations
from bisect import bisect
import random

from structopt.tools import root, single_core, parallel, disjoint_set_merge
from .all_close_atom_positions import all_close_atom_positions
from .diversify_module import diversify_module


class Fingerprinters(object):
    """ """

    @single_core
    def __init__(self, parameters):
        self.parameters = parameters
        self.fingerprinters = {getattr(self, name): self.parameters[name]['probability'] for name in self.parameters}
        self.kwargs = {getattr(self, name): self.parameters[name]['kwargs'] for name in self.parameters}
        total_probability = sum(self.fingerprinters.values())
        self.fingerprinters[None] = 1.0 - total_probability
        self.selected_fingerprinter = None

    @single_core
    def select_fingerprinter(self):
        # Implementation from https://docs.python.org/3/library/random.html -- Ctrl+F "weights"
        choices, weights = zip(*self.fingerprinters.items())
        cumdist = list(accumulate(weights))
        x = random.random() * cumdist[-1]
        self.selected_fingerprinter = choices[bisect(cumdist, x)]

    @single_core
    def remove_duplicates(self, population, nkeep):
        assert nkeep <= len(population)
        if nkeep == len(population):
            return

        kwargs = self.kwargs[self.selected_fingerprinter]
        equivalent_pairs = []
        for pair in combinations(population, 2):
            are_the_same = self.selected_fingerprinter(*pair, **kwargs)
            if are_the_same:
                equivalent_pairs.append(pair)

        if equivalent_pairs:
            ids = [i.id for i in population]
            equivalent_pairs = [(a.id, b.id) for a, b in equivalent_pairs]
            equivalent_individuals = disjoint_set_merge(ids, equivalent_pairs)
            to_keep = []
            for individuals in equivalent_individuals:
                to_keep.append(random.choice(tuple(individuals)))

            while len(to_keep) < nkeep:
                to_keep.append(
                    random.choice(
                        [id for id in ids if id not in to_keep]
                    )
                )

            new_population = [population[id] for id in to_keep]
            population.replace(new_population)

    @single_core
    def post_processing(self):
        pass

    @staticmethod
    @functools.wraps(all_close_atom_positions)
    def all_close_atom_positions(individual1, individual2, **kwargs):
        return all_close_atom_positions(individual1, individual2, **kwargs)

    @staticmethod
    @functools.wraps(diversify_module)
    def diversify_module(individual1, individual2, **kwargs):
        return diversify_module(individual1, individual2, **kwargs)

