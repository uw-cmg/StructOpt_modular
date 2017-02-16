import functools
import random, time
from itertools import accumulate, combinations
from bisect import bisect
from mpi4py import MPI
from structopt.tools import root, single_core, parallel, disjoint_set_merge
from structopt.tools.parallel import allgather
import gparameters
from .all_close_atom_positions import all_close_atom_positions
from .diversify_module import diversify_module


class Fingerprinters(object):
    """ """
    kwargs = ['keep_best']

    @single_core
    def __init__(self, parameters):
        self.parameters = parameters
        self.fingerprinters = {getattr(self, name): self.parameters[name]['probability'] for name in self.parameters if name not in self.kwargs}
        self.function_kwargs = {getattr(self, name): self.parameters[name]['kwargs'] for name in self.parameters if name not in self.kwargs}
        total_probability = sum(self.fingerprinters.values())
        self.fingerprinters[None] = 1.0 - total_probability
        self.selected_fingerprinter = None
        self.equivalent_pairs = []

    @single_core
    def select_fingerprinter(self):
        # Implementation from https://docs.python.org/3/library/random.html -- Ctrl+F "weights"
        choices, weights = zip(*self.fingerprinters.items())
        cumdist = list(accumulate(weights))
        x = random.random() * cumdist[-1]
        self.selected_fingerprinter = choices[bisect(cumdist, x)]

    @single_core
    def remove_duplicates(self, population, nkeep, keep_best=True):
        assert nkeep <= len(population)
        if nkeep == len(population):
            return

        # Apply the selected fingerprinter to all pairs of individuals
        self.get_equivalent_pairs(population)

        if self.equivalent_pairs and gparameters.mpi.rank == 0:
            if keep_best:
                fits = [(individual.fitness, individual.id) for individual in population]
                best = sorted(fits)[0][1]

            ids = [i.id for i in population]
            self.equivalent_pairs = [(a.id, b.id) for a, b in self.equivalent_pairs]
            equivalent_individuals = disjoint_set_merge(ids, self.equivalent_pairs)
            to_keep = set()
            killed = set()
            for equivalent in equivalent_individuals:
                if keep_best and best in equivalent:
                    to_keep.add(best)
                    for x in equivalent:
                        if x is not best:
                            killed.add(x)
                else:
                    equivalent = sorted(equivalent, key=lambda id: population[id].fitness)
                    to_keep.add(equivalent[0])
                    for x in equivalent[1:]:
                        killed.add(x)

            while len(to_keep) < nkeep:
                to_keep.add(
                    random.choice(
                        [id for id in ids if id not in to_keep]
                    )
                )

            killed = [population[id] for id in killed]
            new_population = [population[id] for id in to_keep]
            population.replace(new_population)
            print('Finished killing in {}s'.format(time.time() - t0))
            return killed
        else:
            return []

    @parallel
    def get_equivalent_pairs(self, population):
        """Get equivalent pairs in parallel

        Args:
            population (Population): the population
        """
        kwargs = self.function_kwargs[self.selected_fingerprinter]
        rank = gparameters.mpi.rank

        ncores = gparameters.mpi.ncores
        pairs_per_core = {r: [] for r in range(ncores)}
        for i, pair in enumerate(combinations(population, 2)):
            pairs_per_core[i % ncores].append(pair)
        for pair in pairs_per_core[rank]: 
            are_the_same = self.selected_fingerprinter(*pair, **kwargs)
            if are_the_same:
                self.equivalent_pairs.append(pair)
        print("Finished running {} pairs on rank {}".format(len(pairs_per_core[rank]), rank))
        return None
    
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

