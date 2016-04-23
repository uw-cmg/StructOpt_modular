import functools.wraps
import structopt.globals._mutations


class Mutations(structopt.globals.mutations.Mutations):
    """ """
    @staticmethod
    @functools.wraps(structopt.globals._mutations.add_atoms)
    def add_atoms(individual):
        # Do something interesting to the individual here
        individual = structopt.globals._mutations.add_atoms(individual)
        # Undo that interesting thing to the individual here
        return individual
