import numpy as np

from structopt.tools import root, single_core, parallel

@single_core
def roulette(population, fits, nkeep):
    """Select individuals with a probability proportional to their fitness.
    Fitnesses are renormalized from 0 - 1.
    """

    # Normalize fits from 0 (min fit) to 1 (max fit)
    fit_max = max(fits)
    fits = np.array([-(fit - fit_max) for fit in fits])
    fits /= np.nan_to_num(max(fits))

    # Generate probabilities and pick individuals
    p = np.nan_to_num(fits / np.sum(fits))

    # If less species have nonzero probability than nkeep, select
    # all nonzero probability and select random zero probability
    indexes_nonzero_p = [i for i, p_i in enumerate(p) if p_i != 0]
    indexes_zero_p = [i for i, p_i in enumerate(p) if p_i == 0]
    if len(indexes_nonzero_p) < nkeep:
        n_zero_p_to_add = nkeep - len(indexes_zero_p) - len(indexes_nonzero_p)
        indexes_zero_p_keep = np.random.choice(indexes_zero_p, n_zero_p_to_add, replace=False)
        indexes_keep = np.append(indexes_nonzero_p, indexes_zero_p_keep)
    else:
        indexes_keep = np.random.choice(len(population), nkeep, replace=False, p=p)    

    new_population = [population[i] for i in indexes_keep]
    population.replace(new_population)

    return None
