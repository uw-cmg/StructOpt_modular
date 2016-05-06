import random
from operator import attrgetter
from itertools import accumulate


def cost(population, fits, nkeep):
    """Selection function to order and select structures based on simple relative cost"""
    """The cost selection option is a roulette wheel selection method.  Each individual is given a selection probability relative to the least fit individual in the population.  This method will allow duplicate structures to exist within the selection list but will prevent duplicates from being selected adjacently, which will prevent a structure from attempting a crossover with itself.  The primary difference between the cost selection method and the rank selection method is the weighting each individual receives.  In the cost selection method, the weight is assigned by not only giving those with higher fitness a higher weight, but by giving the higher fitness individuals a proportionally higher weight.  This can be useful in populations with higher diversity, as it will proportionally favor those individuals with higher fitness better than the rank selection method."""
    """Modifies `population` in place."""
    # TODO I don't actually think that this is true: This method will allow duplicate structures to exist within the selection list but will prevent duplicates from being selected adjacently.

    # Sort by fitnesses and calculate the cumulative probability of selecting each individual
    sorted_fits, sorted_population = zip(*sorted(zip(fits, population)))
    population.replace(list(sorted_population))

    try:
        norms = [fit - fits[nkeep+1] for fit in fits[0:nkeep]]
    except IndexError:
        norms = [fit - fits[nkeep-1] for fit in fits[0:nkeep]]
    cumulative_probability = list(accumulate(norms))

    chosen = []
    for i in range(nkeep):
        rand = random.random()
        counter = 0
        while cumprob[counter] < rand:
            counter += 1
        if population[counter] in chosen:
            a = random.choice(population)
            chosen.append(a)
        else:
            chosen.append(population[counter])

    population.replace(chosen)

