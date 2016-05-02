import functools.wraps
import structopt.individual.mutations


class Mutations(structopt.individual.mutations.Mutations):
    """ """
    @staticmethod
    @functools.wraps(structopt.individual.mutations.add_atoms)
    def add_atoms(individual):
        # Do something interesting to the individual here
        individual = structopt.individual.mutations.add_atoms(individual)
        # Undo that interesting thing to the individual here
        return individual

