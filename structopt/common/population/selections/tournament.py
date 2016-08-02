import scipy
import numpy as np
from copy import deepcopy

def tournament(population, fits, tournament_size=5,
               unique_pairs=False, unique_parents=False):
    """Selects pairs in seperate "tournaments", where a subset of the
    population are randomly selected and the highest fitness allowed to pass.
    In addition to a population, their fits, and end population size, takes in
    a tournament size parameter.

    Parameters
    ----------
    population : StructOpt.Population object
        The population of individuals needed to be trimmed
    fits : list
        List of fitnesses that correspond to the population.
    nkeep : int
        The number of individuals in the next popoulation
    tournament_size : int
        The number of individuals in each tournament. If 1,
        tournament is the same as random selection. If
        len(population), corresponds to the "best" selection process
    unique_pairs : bool
        If True, all combinations of parents are unique.
        True increases the diveristy of the population.
    unique_parents : bool
        If True, all parents can only mate with on other individual.
        True increases the diversity of the population.


    Output
    ------
    out : None
        Population modified in place.
    """

    # Get ranks of each population value based on its fitness
    ranks = scipy.stats.rankdata(fits, method='ordinal')
    inds_population = [individual.id for individual in population]
    n_pairs = int(len(fits) / 2)    
    pairs_ind = []

    for i in range(n_pairs):
        # Perform tournament for father selection
        if len(inds_population) > tournament_size:
            tournament_father = np.random.choice(inds_population, tournament_size)
        else:
            tournament_father = inds_population
        max_rank = min([ranks[i] for i in tournament_father])
        ind_father = list(ranks).index(max_rank)

        # Choose the second parent based on renormalized probabilities
        # First remove the father from the population and make a temporary
        # population for choosing the mother. If unique_parents
        # is on, remove him from the global population as well.
        inds_population_temp = deepcopy(inds_population)
        ind_delete = inds_population.index(ind_father)
        del inds_population_temp[ind_delete]

        if unique_parents:
            del inds_population[ind_delete]

        # Now remove mothers that would make repeat father/mother pairs
        if unique_pairs:
            for pair in deepcopy(pairs_ind):
                if ind_father in pair:
                    del pair[pair.index(ind_father)]
                    ind_mother = pair[0]
                    ind_delete = inds_population_temp.index(ind_mother)
                    del inds_population_temp[ind_delete]

        # Using subset, perform tournament to find mother
        if len(inds_population_temp) > tournament_size:
            tournament_mother = np.random.choice(inds_population_temp, tournament_size)
        else:
            tournament_mother = inds_population_temp
        max_rank = min([ranks[i] for i in tournament_mother])
        ind_mother = list(ranks).index(max_rank)

        if unique_parents:
            ind_delete = inds_population.index(ind_mother)
            del inds_population[ind_delete]

        pairs_ind.append([ind_father, ind_mother])

    pairs = [[population[i], population[j]] for i, j in pairs_ind]    

    return pairs

