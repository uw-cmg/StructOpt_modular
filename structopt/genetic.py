import sys, os
import logging

import structopt
from structopt.common.population import Population


class GeneticAlgorithm(object):
    """Defines methods to run a genetic algorithm optimization using the functions in the rest of the library."""

    def __init__(self, population, convergence):
        self.logger = logging.getLogger('default')

        self.population = population
        self.convergence = convergence

        self.generation = 0

        # Set starting convergence
        self.converged = False


    def run(self):
        if logging.parameters.rank == 0:
            print("Starting main Opimizer loop!")
        while not self.converged:
            self.step()
        if logging.parameters.rank == 0:
            print("Finished!")


    def step(self):
        if logging.parameters.rank == 0:
            print("Starting generation {}".format(self.generation))
        sys.stdout.flush()
        if self.generation > 0:
            fits = [individual._fitness for individual in self.population]
            parents = self.population.select(fits)
            children = self.population.crossover(parents)
            self.population.extend(children)
            mutated_population = self.population.mutate()
            self.population.replace(mutated_population)
        self.population.relax()
        fits = self.population.fitness()
        self.population.kill(fits)
        self.check_convergence()
        if logging.parameters.rank == 0:
            self.post_processing_step()
        self.population.generation += 1
        self.generation += 1


    def check_convergence(self):
        if self.generation >= self.convergence.max_generations:
            self.converged = True
        else:
            self.converged = False


    def post_processing_step(self):
        # Save the fitnesses for each individual
        fitness_logger = logging.getLogger('fitness')
        for individual in self.population:
            fitness_logger.info('Generation {}: {}'.format(self.generation, ', '.join([str(x) for x in individual.fits])))

        # Save the XYZ file for each individual
        for individual in self.population:
            path = 'XYZs'
            if not os.path.exists(path):
                os.makedirs(path)
            path = os.path.join(path, 'generation{}'.format(self.generation))
            if not os.path.exists(path):
                os.makedirs(path)
            individual.write(os.path.join(path, 'individual{}.xyz'.format(individual.index)))


    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass


if __name__ == "__main__":
    import sys
    import structopt
    import random
    import numpy as np
    from genetic import GeneticAlgorithm

    parameters = structopt.setup(sys.argv[1])
    random.seed(parameters.seed)
    np.random.seed(parameters.seed)
    
    population = Population(parameters=parameters)

    with GeneticAlgorithm(population=population,
                                    convergence=parameters.convergence
                                    ) as optimizer:
        optimizer.run()

