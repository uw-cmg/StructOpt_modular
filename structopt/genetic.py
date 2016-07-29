import sys, os
import logging
import shutil

import structopt
from structopt.common.population import Population


class GeneticAlgorithm(object):
    """Defines methods to run a genetic algorithm optimization using the functions in the rest of the library."""

    def __init__(self, population, convergence, post_processing):
        self.logger = logging.getLogger('default')

        self.population = population
        self.convergence = convergence
        self.post_processing = post_processing

        self.generation = 0
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
        self.population.sort()
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
            line = 'Generation {}, Individual {}:'.format(self.generation, individual.index)
            for module in individual.fits:
                line += ' {}: {}'.format(module, individual.fits[module])
            fitness_logger.info(line)

        # Save the XYZ file for each individual
        for individual in self.population:
            path = os.path.join(logging.parameters.path, 'XYZs')
            os.makedirs(path, exist_ok=True)
            path = os.path.join(path, 'generation{}'.format(self.generation))
            os.makedirs(path, exist_ok=True)
            individual.write(os.path.join(path, 'individual{}.xyz'.format(individual.index)))

        self.clear_XYZs()

        # Save the genealogy
        tags = ['' for _ in self.population]
        for i, individual in enumerate(self.population):
            tags[i] = '{ctag}{index}{mtag}'.format(ctag=individual.crossover_tag or '', index=individual.index, mtag=individual.mutation_tag or '')
            individual.crossover_tag = None
            individual.mutation_tag = None
        genealogy_logger = logging.getLogger('genealogy')
        genealogy_logger.info(' '.join(tags))

    def clear_XYZs(self):
        """Depending on the value in the post_processing dictionary, clear old
        XYZ files to save space. Specified by the parameters.post_processing.XYZs kwarg.
        Takes an integer n  value. Behavior depends on sign of integer. Always 
        includes the first and last generation

        -n : Only generation up to current generation - n are kept

        n : Every n generation is kept. """

        path = None

        n = self.post_processing['XYZs']
        assert(type(n) is int)

        # If nothing is getting removed
        if self.generation == 0:
            return

        # Keeping the last n generations
        if n < 0 and self.generation > -n:
            path = os.path.join(logging.parameters.path, 'XYZs/generation{}'.format(self.generation + n))
        # Keeping every n generation
        elif n > 0 and self.generation > 1 and self.generation % n != 1:
            path = os.path.join(logging.parameters.path, 'XYZs/generation{}'.format(self.generation - 1))

        if path is not None:
            shutil.rmtree(path)

        return

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
                          convergence=parameters.convergence,
                          post_processing=parameters.post_processing) as optimizer:    
        optimizer.run()

