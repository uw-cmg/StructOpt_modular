from copy import deepcopy
import numpy as np
import scipy

def rank(population, fits, p_min=None, unique_pairs=False, unique_parents=False):
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
    unique_pairs : bool
        If True, all combinations of parents are unique.
        True increases the diveristy of the population.
    unique_parents : bool
        If True, all parents can only mate with on other individual.
        True increases the diversity of the population.

    Returns
    -------
    out : list
        A list of pairs of crossover pairs. Is always at most half the size
        of the population.
    """

    # Get ranks of each population value based on its fitness
    ranks = scipy.stats.rankdata(fits, method='ordinal')

    # Work with indexes of the population instead of the population
    inds_population = list(range(len(population)))

    # Get probabilities based on linear ranking
    if p_min is None:
        p_min = 1.0 / len(population) ** 2
    N = len(fits)
    eta_min = p_min * N
    eta_max = 2 - eta_min
    p_max = eta_max / N
    p = p_min + (p_max - p_min)*(N - ranks)/(N - 1)

    # Construct list of parents
    n_pairs = int(len(fits) / 2)
    pairs_ind = []

    for i in range(n_pairs):
        # Choose the first parent based on probabilties
        ind_father = np.random.choice(inds_population, p=p)

        # Choose the second parent based on renormalized probabilities
        # First remove the father from the population and make a temporary
        # population for choosing the mother. If unique_parents
        # is on, remove him from the global population as well.
        inds_population_temp, p_temp = list(inds_population), list(p)
        ind_delete = inds_population.index(ind_father)
        del inds_population_temp[ind_delete]
        p_temp = np.delete(p_temp, ind_delete)

        if unique_parents:
            del inds_population[ind_delete]
            p = np.delete(p, ind_delete)
            p /= sum(p)

        # Now remove mothers that would make repeat father/mother pairs
        if unique_pairs:
            for pair in deepcopy(pairs_ind):
                if ind_father in pair:
                    del pair[pair.index(ind_father)]
                    ind_mother = pair[0]
                    ind_delete = inds_population_temp.index(ind_mother)
                    del inds_population_temp[ind_delete]
                    p_temp = np.delete(p_temp, ind_delete)

        p_temp /= sum(p_temp)
        ind_mother = np.random.choice(inds_population_temp, p=p_temp)

        if unique_parents:
            ind_delete = inds_population.index(ind_mother)
            del inds_population[ind_delete]
            p = np.delete(p, ind_delete)
            p /= sum(p)

        pairs_ind.append([ind_father, ind_mother])

    # Construct the parents from the indices
    pairs = [[population[i], population[j]] for i, j in pairs_ind]

    return pairs
