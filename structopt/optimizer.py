class Optimizer(object):
    def __init__(self):
        # Initialize random number seed
        random.seed(self.seed)

        # Prep output monitoring

        # Set starting convergence
        self.converged = False

    def run(self):
        while not self.converged:
            self.step()

    def step(self):
        pass
        self.check_convergence()

    def check_convergence(self):
        self.converged = False

