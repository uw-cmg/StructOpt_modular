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
    c_indiv = individual.positions.sum(0)/natoms
    c_swarm = best_swarm.positions.sum(0)/natoms
    c_particle = best_particle.positions.sum(0)/natoms

    delta_particle = numpy.zeros([natoms, 3])
    delta_swarm = numpy.zeros([natoms, 3])
    
    swarm = best_swarm.positions - c_swarm + c_indiv
    particle = best_particle.positions - c_particle + c_indiv
    
    for i in range(natoms):
        min_dist = 1000
        jj = 0
        for j in range(len(swarm)):
            dist = numpy.linalg.norm(swarm[j] - individual.positions[i])
            if dist < min_dist:
                delta_swarm[i] = swarm[j] - individual.positions[i]
                min_dist = dist
                jj = j
        numpy.delete(swarm, jj)

    for i in range(natoms):
        min_dist = 1000
        jj = 0
        for j in range(len(particle)):
            dist = numpy.linalg.norm(particle[j] - individual.positions[i])
            if dist < min_dist:
                delta_particle[i] = particle[j] - individual.positions[i]
                min_dist = dist
                jj = j
        numpy.delete(particle, jj)    

    individual.set_velocities(omega * individual.velocities +
                              phi_p * rp * delta_particle  +
                              phi_g * rg * delta_swarm)
    individual.set_positions(individual.positions + individual.velocities)
    
    return None

