import functools

import structopt
from structopt.tools import root, single_core, parallel


class Pso_moves(object):
    """ """

    @single_core
    def __init__(self, parameters):
        self.parameters = parameters


    @single_core
    def move(self, population, best_swarm, best_particles):
        for i, individual in enumerate(population):
            individual.pso_moves.update_particle(individual, best_swarm, best_particles[i])
        return population

    @single_core
    def post_processing(self):
        pass

