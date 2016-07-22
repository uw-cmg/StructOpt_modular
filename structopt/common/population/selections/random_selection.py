import random as random
from itertools import combinations

def random_selection(population, fits):
    pairs = [pair for pair in combinations(population, 2)]
    return pairs

