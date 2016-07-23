import random as random
from itertools import combinations

def random_selection(population, fits, prob):
    pairs = []
    for pair in combinations(population, 2):
        if random.random() < prob:
            pairs.append(pair)
    return pairs

