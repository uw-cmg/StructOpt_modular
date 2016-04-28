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
        population.crossover()
        for individual in self.population:
            individual.mutate()
        for individual in self.population:
            individual.relax()
        for individual in self.population:
            individual.fitness()
        population.kill()
        population.select()
        self.check_convergence()

    def check_convergence(self):
        self.converged = False

