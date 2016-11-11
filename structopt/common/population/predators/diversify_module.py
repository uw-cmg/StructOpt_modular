import random


def diversify_module(population, fits, nkeep, module, min_diff):
    """This predator eliminates individuals with similar module values"""

    module_fits = [getattr(individual, module) for individual in population]
    ids = [individual.id for individual in population]
    module_fits, ids = zip(*sorted(list(zip(module_fits, ids))))
    diffs = [j - i for i, j in zip(module_fits[:-1], module_fits[1:])]
    ids_delete = []
    best_id = ids[list(fits).index(min(fits))]
    for i in range(len(diffs)):
        if diffs[i] < min_diff:
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
        population.replace(list(sorted_population[:nkeep]))
