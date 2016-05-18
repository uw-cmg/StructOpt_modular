
def best(population, fits, nkeep):
    # Sort by fitnesses and calculate the cumulative probability of selecting each individual
    sorted_fits, sorted_population = zip(*sorted(zip(fits, population), key=lambda pair: pair[0]))
    population.replace(list(sorted_population[:nkeep]))

