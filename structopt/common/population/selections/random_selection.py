from itertools import combinations
import random


def random_selection(population, fits):
    pairs = [pair for pair in combinations(population, 2)]
    random.shuffle(pairs)
    pairs = pairs[:len(population) // 2]
    return pairs

