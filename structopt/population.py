import structopt


class Population(list):
    def __init__(self):
        structopt.parameters.global.structure_type

        # Generate/load initial structures
        pass

        # Allow structopt to see the population
        structopt.population = self

    def crossover(self):
        self.crossovers.select_crossover()
        return self.crossovers.crossover(self.individuals)

    def select(self):
        self.selections.select_selection()
        return self.selections.select(self.individuals)

    def kill(self):
        self.predators.select_predator()
        return self.predators.kill(self.individuals)

