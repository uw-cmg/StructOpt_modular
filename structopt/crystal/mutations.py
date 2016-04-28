import functools.wraps
import structopt.globals.individual._mutations
import structopt.globals.individual.mutations


class Mutations(structopt.globals.individual.mutations.Mutations):
    """ """
    @staticmethod
    @functools.wraps(structopt.globals.individual._mutations.add_atoms)
    def add_atoms(individual):
        # Do something interesting to the individual here
        individual = structopt.globals.individual._mutations.add_atoms(individual)
        # Undo that interesting thing to the individual here
        return individual
