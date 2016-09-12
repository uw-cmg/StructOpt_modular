import random

def permutation(individual):
    """Move function to perform Permutation of one atom based on atomlist
    Inputs:
        indiv = Individual class object to be altered
        Optimizer = Optimizer class object with needed parameters
    Outputs:
        indiv = Altered Individual class object
    """

    a1 = individual[random.randint(0, individual.get_number_of_atoms() - 1)]
    opts = [i for i in individual if i.symbol != a1.symbol]
    a2 = opts[random.randint(0,len(opts) - 1)]
    a1.symbol, a2.symbol = a2.symbol, a1.symbol

    return
