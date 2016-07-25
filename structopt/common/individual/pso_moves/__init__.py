import logging
import functools
from collections import defaultdict

import structopt
from .update_particle import update_particle
from structopt.tools import root, single_core, parallel


class Pso_moves(object):
    """ """

    @single_core
    def __init__(self, parameters):
        # These variables never change
        self.parameters = parameters
        self.kwargs = defaultdict(dict)
        self.kwargs.update( {getattr(self, name): kwords for name, kwords in self.parameters.kwargs.items()} )


    @single_core
    def move(self, individual, best_swarm, best_particle):
        logger = logging.getLogger("default")        
        individual._relaxed = False
        individual._fitted = False
        return self.update_particle(individual, best_swarm, best_particle)


    @single_core
    def post_processing(self):
        pass


    @staticmethod
    @functools.wraps(update_particle)
    def update_particle(individual, best_swarm, best_particle):
        return update_particle(individual,
                               best_swarm,
                               best_particle,
                               omega=self.parameters.omega,
                               phi_p=self.parameters.phi_p,
                               phi_g=self.parameters.phi_g)

