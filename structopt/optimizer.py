import sys
import random
import logging

import structopt
from structopt.common.population import Population


class Optimizer(object):
    __version__  = 'StructOpt_v3.0'

    def __init__(self, parameters):
        self.logger = logging.getLogger('default')

        # Get parameters from StructOpt space
        self.parameters = parameters

        # Initialize random number seed
        random.seed(self.parameters.seed)

        self.generation = 0

        # Create the population (not ready yet)
        self.population = Population(parameters=self.parameters)

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
        self.generation += 1


    def check_convergence(self):
        if self.generation >= self.parameters.convergence.maxgen:
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

    parameters = structopt.setup(sys.argv[1])

    with structopt.Optimizer(parameters) as optimizer:
        optimizer.run()

