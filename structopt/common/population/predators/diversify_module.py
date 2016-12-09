import random


def diversify_module(population, fits, nkeep, module='LAMMPS', min_diff=0.0001):
    """This module eliminates individuals that have similar fitness values.
    The goal is to diversify the population and eliminate duplicate 
    individuals. It is often useful to include this in the optimizer along
    with another predator.

    Parameters
    ----------
    population : Population
        A population of individuals
    fits : list
        Fitnesses that corresponds to population
    nkeep : int
        The number of individuals to keep. In a GA run, corresponds
        to the sum of each generators number_of_individuals
    module : str
        The name of the module for evaluating similar fitnesses.
    min_diff : float
        The minimum difference between individual fitness values for judging
        they are the same.
    """

    module_fits = [getattr(individual, module) for individual in population]
    ids = [individual.id for individual in population]
    best_id = ids[list(fits).index(min(fits))]
    module_fits, ids = zip(*sorted(list(zip(module_fits, ids))))
    diffs = [j - i for i, j in zip(module_fits[:-1], module_fits[1:])]
    ids_delete = []

    for i, diff in enumerate(diffs):
        if diff < min_diff:
            ids_delete.append(ids[i + 1])

    if best_id in ids_delete:
        del ids_delete[ids_delete.index(best_id)]

    if len(population) - len(ids_delete) <= nkeep:
        random.shuffle(ids_delete)
        ids_delete = ids_delete[:len(population) - nkeep]
        new_population = [individual for individual in population if individual.id not in ids_delete]
        population.replace(new_population)
    else:
        new_population = [individual for individual in population if individual.id not in ids_delete]
        new_module_fits = [individual._fitness for individual in new_population]
        sorted_module_fits, sorted_population = zip(*sorted(zip(new_module_fits, new_population), key=lambda pair: pair[0]))
        sorted_population = list(sorted_population)
        new_population = [sorted_population.pop(0)]
        random.shuffle(sorted_population)
        new_population += sorted_population[:nkeep-1]
        population.replace(new_population)
