import sys
import random
import logging

import structopt
from structopt.common.population import Population


def update_particle(individual, best_swarm, best_particle, omega, phi_p, phi_g):
    natoms = len(individual.positions)
    
    individual.velocities = omega * individual.velocities
                            + phi_p*rp*(best_particle.positions - individual.positions) 
                            + phi_g*rg*(best_swarm.positions - individual.positions)
    individual.positions += individual.velocities

    return None

