import random

import structopt
from structopt.common.population import Population


class Optimizer(object):
    def __init__(self):
        # Initialize random number seed
        random.seed(structopt.parameters.globals.seed)

        self.generation = 0

        # Create the population
        self.population = Population()

        # Prep output monitoring

        # Set starting convergence
        self.converged = False

    def run(self):
        while not self.converged:
            self.step()

    def step(self):
        child1, child2 = self.population.crossover()
        if child1 is not None:
            self.append(child1)
        if child2 is not None:
            self.append(child2)
        self.population.mutate()
        self.population.relax()
        fits = self.population.fitness()
        self.population.kill()
        self.population.select(fits, structopt.parameters.globals.total_number_of_individuals)
        self.check_convergence()
        self.generation += 1

    def check_convergence(self):
        if self.generation > 10:
            self.converged = True
        else:
            self.converged = False


if __name__ == "__main__":
    import sys
    import structopt

    structopt.setup(sys.argv[1])

    optimizer = structopt.Optimizer()
    optimizer.run()
