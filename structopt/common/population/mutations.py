import numpy as np

from structopt.tools import root, single_core, parallel


class Mutations(object):
    """ """

    @single_core
    def __init__(self, parameters):
        self.parameters = parameters
        if 'preserve_best' in parameters and parameters['preserve_best'] == True:
            self.preserve_best = True
        else:
            self.preserve_best = False

    @single_core
    def mutate(self, population):
        fits = [getattr(ind, '_fitness') if ind._fitted else np.inf for ind in population]
        min_fit_index = fits.index(np.amin(fits))
        for i, individual in enumerate(population):
            individual.mutations.select_mutation()
            if self.preserve_best and i == min_fit_index:
                individual.mutations.selected_mutation = None
            if individual.mutations.selected_mutation is not None:
                individual.mutate()
        return population


    @single_core
    def post_processing(self):
        pass
