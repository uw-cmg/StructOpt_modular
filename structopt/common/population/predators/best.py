
def best(population, fits, nkeep):
    """Sorts individuals by fitness and keeps the top nkeep fitnesses.
    
    Parameters
    ----------
    population : Population
        A population of individuals
    fits : list
        Fitnesses that corresponds to population
    nkeep : int
        The number of individuals to keep. In a GA run, corresponds
        to the sum of each generators number_of_individuals
    """
    sorted_fits, sorted_population = zip(*sorted(zip(fits, population), key=lambda pair: pair[0]))
    population.replace(list(sorted_population[:nkeep]))

