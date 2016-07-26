'''This contains a StructOpt calculator for submitting jobs to queue, tracking their progress, and reading their output. Currently only works with genetic.py'''

import os
import re

import numpy as np

class StructOpt(object):

    def __init__(self, calcdir=None, optimizer=None, parameters=None):

        if calcdir is None:
            self.calcdir = os.getcwd()
        else:
            self.calcdir = os.path.expandvars(calcdir)
        self.cwd = os.getcwd()
        self.optimizer = optimizer
        self.parameters = parameters

        self.fitness = None

        self.system_name = os.path.basename(self.calcdir)

    def __enter__(self):
        '''On enter, make sure directory exists. Create it if necessary
        and change into the directory. Then return the calculator.'''

        # Make directory if it doesn't already exist
        if not os.path.isdir(self.calcdir):
            os.makedirs(self.calcdir)

        # Now change into the new working directory
        os.chdir(self.calcdir)
        self.initialize()

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        '''On exit, change back to the original directory'''

        os.chdir(self.cwd)

        return

    def initialize(self):
        '''Gets the status of the calculation. Meant to be run when within
        the calculation directory'''

        pass

    def read_fitness(self):
        '''Reads fitness.log and stores the data'''

        all_fitness = []
        pattern = '.* Generation (.*), Individual (.*): (.*)'

        current_population = []
        current_generation = 0
        with open('fitnesses.log') as fitness_file:
            for line in fitness_file:
                match = re.match(pattern, line, re.I|re.M)
                if match:
                    generation = int(match.group(1))
                    individual = int(match.group(2))
                    fitness = float(match.group(3))
                else:
                    continue

                if generation > current_generation:
                    current_generation = generation
                    all_fitness.append(current_population)
                    current_population = []

                if len(current_population) < individual + 1:
                    current_population += [None]*(individual + 1 - len(current_population))

                current_population[individual] = fitness

        self.fitness = np.array(all_fitness)

        return

    def get_fitnesses(self):
        """Returns a list of fitnesses of all individuals in all generations.

        Output
        ------
        out : N x M list
            N = Number of generations
            M = Number of individuals in each generation
            out[i][j] gives fitness of individual j in generation i.
        """

        if self.fitness is None:
            self.read_fitness()

        return self.fitness

    def get_avg_fitnesses(self):
        """Returns a list of the average fitness of each generation

        Output
        ------
        out : N list
            N = Number of generations
            out[i] gives average fitness of generation i
        """

        if self.fitness is None:
            self.read_fitness()

        return np.average(self.fitness, axis=1)
