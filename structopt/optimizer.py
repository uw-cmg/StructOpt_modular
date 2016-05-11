import random

import structopt
from structopt.common.population import Population


class Optimizer(object):
    def __init__(self):
        
        ### self.logger = logging.getLogger('default')

        # Get parameters from StructOpt space
        setattr(self, 'parameters', structopt.parameters)
        
        # Initialize random number seed
        random.seed(self.parameters.globals.seed)

        self.generation = 0

        # Create the population (not ready yet)
        ### self.population = Population() 
    

        # Prep output monitoring

        # Set starting convergence
        self.converged = False

    def run(self):
        while not self.converged:
            self.step()

    def step(self):
        self.population.crossover()
        self.population.mutate()
        self.population.relax()
        fits = self.population.fitness()
        self.population.kill()
        self.population.select(fits)
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
