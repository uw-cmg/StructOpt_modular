import sys
import random
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

        # Prep output monitoring

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
            pairs = self.population.select(fits)
            crossed_population = self.population.crossover(pairs)
            self.population.replace(crossed_population)
            mutated_population = self.population.mutate()
            self.population.replace(mutated_population)
        self.population.relax()
        fits = self.population.fitness()
        self.population.kill(fits)
        self.check_convergence()
        self.population.generation += 1
        self.generation += 1


    def check_convergence(self):
        if self.generation >= self.convergence.maxgen:
            self.converged = True
        else:
            self.converged = False
    

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass


if __name__ == "__main__":
    import sys
    import structopt
    import numpy as np

    parameters = structopt.setup(sys.argv[1])

    random.seed(parameters.seed)
    np.random.seed(parameters.seed)

    population = Population(parameters=parameters)

    with structopt.GeneticAlgorithm(population=population,
                                    convergence=parameters.convergence
                                    ) as optimizer:
        optimizer.run()

