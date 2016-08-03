import random
import scipy
import numpy as np

def tournament(population, fits, nkeep, tournament_size=5):
    """Selects individuals in seperate "tournaments", where a subset of the
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

    Output
    ------
    out : None
        Population modified in place.
    """

    # Get ranks determined by their fitness
    ranks = list(scipy.stats.rankdata(fits, method='ordinal'))

    # Work with indices of the indices of the population
    ids_population = [individual.id for individual in population]
    new_population = []

    fits_keep, ids_keep = [], []
    
    # Run nkeep tournaments
    for _ in range(nkeep):
        if len(ids_population) > tournament_size:
            tournament_ids = np.random.choice(ids_population,
                                              tournament_size, replace=False)
        else:
            tournament_ids = ids_population
        tournament_ranks = [ranks[ids_population.index(i)] for i in tournament_ids]
        best_rank = min(tournament_ranks)
        tournament_winner = tournament_ids[tournament_ranks.index(best_rank)]
        new_population.append(population[tournament_winner])
        fits_keep.append(best_rank)
        ids_keep.append(tournament_winner)
        del ranks[ids_population.index(tournament_winner)]
        del ids_population[ids_population.index(tournament_winner)]

    population.replace(new_population)

    return [ids_keep, fits_keep]

