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
    ranks = scipy.stats.rankdata(fits, method='ordinal')

    # Work with indices of the indices of the population
    indices_population = list(range(len(population)))
    new_population = []

    # Run nkeep tournaments
    for _ in range(nkeep):
        if len(indices_population) > tournament_size:
            tournament_indices = np.random.choice(indices_population,
                                                  tournament_size, replace=False)
        else:
            tournament_indices = indices_population
        tournament_ranks = [ranks[i] for i in tournament_indices]
        max_fit = min(tournament_ranks)
        tournament_winner = tournament_indices[tournament_ranks.index(max_fit)]
        new_population.append(population[tournament_winner])
        indices_population.pop(indices_population.index(tournament_winner))

    population.replace(new_population)

    return
