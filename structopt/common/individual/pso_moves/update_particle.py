import sys
import random
import logging
import numpy
import structopt
from structopt.common.population import Population


def update_particle(individual, best_swarm, best_particle, omega, phi_p, phi_g):
    natoms = len(individual.positions)
    rp = numpy.random.rand(natoms,3)
    rg = numpy.random.rand(natoms,3)
    individual.set_velocities(omega * individual.velocities +
                            phi_p*rp*(best_particle.positions - individual.positions)  +
                            phi_g*rg*(best_swarm.positions - individual.positions))
    individual.set_positions(individual.positions + individual.velocities)

    return None

