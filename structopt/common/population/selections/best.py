import scipy.stats

def best(population, fits):
    """Deterministic selection function that chooses adjacently
    ranked individuals as pairs."""

    # Get ranks of each population value based on its fitness
    ranks = list(scipy.stats.rankdata(fits, method='ordinal'))
    pairs = []

    for i in range(len(population) // 2):
        rank_father = 2*i + 1
        rank_mother = 2*i + 2
        father = population[ranks.index(rank_father)]
        mother = population[ranks.index(rank_mother)]
        pairs.append([father, mother])

    return(pairs)


