import numpy as np
import scipy


def rank(population, fits, nkeep, p_min=None):
    """Selection function that chooses pairs of structures
    based on linear ranking.

    See "Grefenstette and Baker 1989 Whitley 1989".

    Parameters
    ----------
    population : StructOpt population object
        An object inherited from list that contains
        StructOpt individual objects.
    fits : list
        A list of fitnesses of the population
    p_min : float
        The probability of choosing the lowest ranked individual.
        Given population of size N, this should be below 1/nindiv.
        The probability of selecting rank N (worst) to rank 1 (best)
        increases from p_min to (2/N - p_min) in even, (1/N - p_min)
        increments. Defaults to (1/N)^2.

    Returns
    -------
    out : list
        A list of pairs of crossover pairs. Is always at most half the size
        of the population.
    """

    # Get ranks of each population value based on its fitness
    ranks = scipy.stats.rankdata(fits, method='ordinal')
    ids = [individual.id for individual in population]

    # Get probabilities based on linear ranking
    if p_min is None:
        p_min = 1.0 / len(population) ** 2
    N = len(fits)
    eta_min = p_min * N
    eta_max = 2 - eta_min
    p_max = eta_max / N
    p = p_min + (p_max - p_min)*(N - ranks)/(N - 1)

    indexes_keep = np.random.choice(ids, nkeep, replace=False, p=p)
    new_population = [population[i] for i in indexes_keep]
    population.replace(new_population)

    return

