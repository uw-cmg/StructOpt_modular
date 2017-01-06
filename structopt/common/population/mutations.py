import numpy as np

from structopt.tools import root, single_core, parallel


class Mutations(object):
    """ """

    @single_core
    def __init__(self, parameters):
        self.parameters = parameters
        self.preserve_best = False
        self.keep_original = False
        self.keep_original_best = False
        if 'preserve_best' in parameters and parameters['preserve_best']:
            self.preserve_best = True
        if 'keep_original' in parameters and parameters['keep_original']:
            self.keep_original = True
        if 'keep_original_best' in parameters and parameters['keep_original_best']:
            self.keep_original_best = True

    @single_core
    def mutate(self, population):
        fits = [ind._fitness if ind._fitted else np.inf for ind in population]
        min_fit_index = fits.index(min(fits))
        original_ids = [individual.id for individual in population]
        for i, individual in enumerate(population):
            # Skip individuals that were added in previous mutations and
            # have already been mutated
            if individual.id not in original_ids:
                continue

            individual.mutations.select_mutation()

            if self.preserve_best and i == min_fit_index:
                individual.mutations.selected_mutation = None

            if ((self.keep_original and individual.mutations.selected_mutation is not None)
                or (self.keep_original_best and i == min_fit_index)):

                # Save and delete the original individual's mutation
                selected_mutation = individual.mutations.selected_mutation
                individual.mutations.selected_mutation = None

                # Make a new individual, attach the mutation, add to popualtion
                individual = individual.copy()
                individual.mutations.selected_mutation = selected_mutation
                population.add(individual)

            if individual.mutations.selected_mutation is not None:
                individual.mutate()

        return population


    @single_core
    def post_processing(self):
        pass
