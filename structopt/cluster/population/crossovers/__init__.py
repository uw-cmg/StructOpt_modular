import functools
import structopt.common.population.crossovers

from .rotate import rotate

class Crossovers(structopt.common.population.crossovers.Crossovers):

    @staticmethod
    @functools.wraps(rotate)
    def rotate(individual, max_natoms=0.2, move_CN=8, surf_CN=11):
        return rotate(individual, max_natoms, move_CN, surf_CN)
