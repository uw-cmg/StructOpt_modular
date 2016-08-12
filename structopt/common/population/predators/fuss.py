import random
import numpy as np

def fuss(population, fits, nkeep, nbest=1, fusslimit=10):
    """Fixed uniform selection scheme. Aimed at maintaining diversity
    in the population. In the case where low fit is the highest
    fitness, selects a fitness between min(fits) and min(fits) + fusslimit,
    if the difference between the min(fit) and max(fit) is larger than fusslimit.

    Parameters
    ----------
    population : structopt.Population object
        The population to be shortened
    fits : list
        A list of the fitnesses
    nkeep : int
        The number of individuals to keep
    nbest : int
        The top n individuals to always keep
    fusslimit : float
        Individuals that have fitness fusslimit 
        worse than the max fitness will not be considered
    """

    ids = [individual.id for individual in population]

    # Find min and max fitness
    minf = min(fits)
    maxf = max(fits)
    if abs(maxf-minf) > fusslimit:
            maxf = minf + fusslimit

    # Select random point on fitness line
    pt = random.uniform(minf,maxf)

    # Calculate the distance of each individual's fitness from that point
    distances = np.absolute(np.array(fits) - pt)

    # Sort distances from min to max
    distances_ids = list(zip(distances, ids))
    distances_ids = sorted(distances_ids, key=lambda i: i[0])

    # Select individuals with lowest distance
    ids_keep = list(list(zip(*distances_ids))[1][:nkeep])

    # Always keep the top nbest individuals
    fits_ids = list(zip(fits, ids))
    fits_ids = sorted(fits_ids, key=lambda i: i[0])
    best_ids = []
    for fit, id in fits_ids[:nbest]:
        if id not in ids_keep:
            best_ids.append(id)
            ids_keep = list(np.random.choice(ids_keep, size=len(ids_keep) - 1))

    ids_keep.extend(best_ids)

    new_population = [population[id] for id in ids_keep]
    population.replace(new_population)
    
    return
