import numpy as np
import random
from copy import deepcopy

def roulette(population, fits, unique_pairs=False, unique_parents=False):
    """Selection function that chooses pairs of structures 
    based on their fitness. Fitnesses are normalized from 0 to 1.

    See "Grefenstette and Baker 1989 Whitley 1989".
    
    Parameters
    ----------
    population : StructOpt population object
        An object inherited from list that contains 
        StructOpt individual objects.
    fits : list
        A list of fitnesses of the population
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

    # Work with indexes of the population instead of the population
    inds_population = list(range(len(population)))
    
    # Normalize fits from 0 (min fit) to 1 (max fit)
    fit_max = max(fits)
    fits = np.array([-(fit - fit_max) for fit in fits])
    fits /= np.nan_to_num(max(fits))

    # Generate probabilities and pick individuals
    p = np.nan_to_num(fits / np.sum(fits))

    # Construct list of parents
    n_pairs = int(len(fits) / 2)
    pairs_ind = []

    for i in range(n_pairs):
        # Choose the first parent based on probabilties
        if np.count_nonzero(p) > 0:
            ind_father = np.random.choice(inds_population, p=p)
        else:
            ind_father = np.random.choice(inds_population)

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
            p = np.nan_to_num(p / np.sum(p))

        # Now remove mothers that would make repeat father/mother pairs
        if unique_pairs:
            for pair in deepcopy(pairs_ind):
                if ind_father in pair:
                    del pair[pair.index(ind_father)]
                    ind_mother = pair[0]
                    ind_delete = inds_population_temp.index(ind_mother)
                    del inds_population_temp[ind_delete]
                    p_temp = np.delete(p_temp, ind_delete)

        p_temp = np.nan_to_num(p_temp / np.sum(p_temp))
        if np.count_nonzero(p_temp) > 0:
            ind_mother = np.random.choice(inds_population_temp, p=p_temp)
        else:
            ind_mother = np.random.choice(inds_population_temp)

        if unique_parents:
            ind_delete = inds_population.index(ind_mother)
            del inds_population[ind_delete]
            p = np.delete(p, ind_delete)
            p = np.nan_to_num(p / np.sum(p))

        pairs_ind.append([ind_father, ind_mother])

    # Construct the parents from the indices
    pairs = [[population[i], population[j]] for i, j in pairs_ind]

    return pairs
