import scipy
import numpy as np


def tournament(population, fits, nkeep, tournament_size=5, keep_best=True):
    """Selects individuals in seperate "tournaments", where a subset of the
    population are randomly selected and the highest fitness allowed to pass.
    In addition to a population, their fits, and end population size, takes in
    a tournament size parameter.

    Parameters
    ----------
    population : Population
        The population of individuals needed to be trimmed
    fits : list
        List of fitnesses that correspond to the population.
    nkeep : int
        The number of individuals to keep. In a GA run, corresponds
        to the sum of each generators number_of_individuals
    tournament_size : int
        The number of individuals in each tournament. If 1,
        tournament is the same as random selection. If
        len(population), corresponds to the "best" selection process
    keep_best : bool
        If set to True, the best individual is always included in the
        following generation.
    """

    # Get ranks determined by their fitness
    ranks = list(scipy.stats.rankdata(fits, method='ordinal'))

    # Work with indices of the indices of the population
    ids_population = [individual.id for individual in population]
    new_population = []

    fits_keep, ids_keep = [], []

    if keep_best:
        best_rank = min(ranks)
        best_id = ids_population[ranks.index(best_rank)]
        fits_keep.append(best_rank)
        ids_keep.append(best_id)
        new_population.append(population[best_id])
        del ranks[ids_population.index(best_id)]
        del ids_population[ids_population.index(best_id)]
        nkeep -= 1 

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

