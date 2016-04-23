import functools.wraps
import structopt.globals._mutations


class Mutations(object):
    """ """
    def __init__(self, parameters):
        self.parameters = parameters
        self.mutations = [getattr(self, name) for name in self.parameters.mutation_options]
        self.selected_mutation = None

    def select_mutation(self):
        self.selected_mutation = random.choice(self.mutations)

    def mutate(self, individual):
        return self.selected_mutation(individual)

    @staticmethod
    @functools.wraps(structopt.globals._mutations.add_atoms)
    def add_atoms(individual):
        return structopt.globals._mutations.add_atoms(individual)

    @staticmethod
    @functools.wraps(structopt.globals._mutations.remove_atoms)
    def remove_atoms(individual):
        return structopt.globals._mutations.remove_atoms(individual)

    @staticmethod
    @functools.wraps(structopt.globals._mutations.remove_surface_atoms)
    def remove_surface_atoms(individual):
        return structopt.globals._mutations.remove_surface_atoms(individual)
