class Optimizer(object):
    def __init__(self):
        # Initialize random number seed
        random.seed(self.seed)

        # Create the population
        self.population = Population()

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
        self.population.fitness()
        self.population.kill()
        self.population.select()
        self.check_convergence()

    def check_convergence(self):
        self.converged = False

