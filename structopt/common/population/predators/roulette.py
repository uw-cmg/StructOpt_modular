import random
from itertools import accumulate

from structopt.tools import root, single_core, parallel


@single_core
def roulette(population, fits, nkeep):
    """Selection function to order and select structures based on simple relative cost. Does not allow duplicates."""
    """The cost selection option is a roulette wheel selection method.  Each individual is given a selection probability relative to the least fit individual in the population.  The primary difference between the cost selection method and the rank selection method is the weighting each individual receives.  In the cost selection method, the weight is assigned by not only giving those with higher fitness a higher weight, but by giving the higher fitness individuals a proportionally higher weight.  This can be useful in populations with higher diversity, as it will proportionally favor those individuals with higher fitness better than the rank selection method. Note that if an individual is chosen twice, a random one is chosen instead (I should probably fix that)."""
    """Modifies `population` in place."""
    # TODO Read the end of the doc string; merge docstrings

    # Sort by fitnesses and calculate the cumulative probability of selecting each individual
    sorted_fits, sorted_population = zip(*sorted(zip(fits, population), key=lambda pair: pair[0]))
    population.replace(list(sorted_population))

    try:
        norms = [fit - fits[nkeep+1] for fit in fits[0:nkeep]]
    except IndexError:
        norms = [fit - fits[nkeep-1] for fit in fits[0:nkeep]]
    sumn = sum(norms)
    # TODO one thing I don't like about this is that the when nkeep == len(population), norms[-1] == 0. That means cumprob[-1] == cumprob[-2] and the last individual can never be selected unless it's random chance
    #print(fits)  # TODO
    #print(norms)  # TODO
    if sumn == 0.0:
        return None  # TODO DELETE; this occurs when all the individuals in the population are identical
        raise ZeroDivisionError
    probability = [x/sumn for x in norms]
    cumulative_probability = list(accumulate(probability))
    #print(probability)  # TODO
    #print(cumulative_probability)  # TODO

    chosen = []
    for i in range(nkeep):
        rand = random.random()
        counter = 0
        while cumulative_probability[counter] < rand:
            counter += 1
        if population[counter] in chosen:
            a = random.choice(population)
            while a in chosen:
                a = random.choice(population)
            chosen.append(a)
        else:
            chosen.append(population[counter])

    for i, ind in enumerate(chosen):
        ind.index = i
    population.replace(chosen)

