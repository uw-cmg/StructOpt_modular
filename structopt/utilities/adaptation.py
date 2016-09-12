import numpy as np
import logging
from structopt.tools.dictionaryobject import DictionaryObject

def adapt(adaptations, population, generation):
    """Function that modifies the Population parameters as is that depends on the status
    of the calculation. Currently supports the following dependents

    generation: adapt after a certain generation

    fitness: adapt after a certain average fitness
    """

    # Perform each adaptation
    adaptations_to_del = []
    for i, adaptation in enumerate(adaptations):

        # Check to see if it needs to be done
        if 'condition' not in adaptation:
            continue

        if 'generation' in adaptation['condition']:
            if adaptation['condition']['generation'] != generation:
                continue

        if 'fitness' in adaptation['condition']:
            current_avg_fitness = np.average([individual._fitness for individual in population])
            if not current_avg_fitness < adaptation['condition']['fitness']:
                continue

        # If we get here, perform the adaptation for both
        # the population and each of its individuals
        for move in ['mutations', 'fitnesses', 'relaxations']:
            if move not in adaptation:
                continue

            new_parameters = adaptation[move]
            for module in new_parameters:
                if type(new_parameters[module]) in [dict, DictionaryObject]:
                    new_parameters[module].setdefault('kwargs', {})

            population.parameters.update(new_parameters)
            for individual in population:
                # Note individual.{}_parameters is not plural
                setattr(individual, '{}_parameters'.format(move[:-1]), new_parameters)
                individual.load_modules()

            print('Setting adaptation: ', new_parameters)

    # Finally delete the adaptations that have been set
    for adaptation in reversed(sorted(adaptations_to_del)):
        del adaptations[adaptation]

    return
