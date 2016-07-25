import sys
import random
import logging

import structopt
from structopt.common.population import Population


def update_particle(individual, best_swarm, best_particle, omega, phi_p, phi_g):
    natoms = len(individual.positions)
    #TODO should use numpy array manipulation instead of for loops 
    for n in range(natoms):
        for j in range(3):
            rp, rg = random.random(), random.random()
            #omega = random.random()
            individual.velocities[n][j] = omega * individual.velocities[n][j] + \
            phi_p*rp*(best_particle.positions[n][j] - individual.positions[n][j]) + \
            phi_g*rg*(best_swarm.positions[n][j] - individual.positions[n][j])

            individual.positions[n][j] += individual.velocities[n][j]

    return None

