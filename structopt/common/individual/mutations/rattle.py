import numpy as np
from structopt.common.crossmodule import get_avg_radii

def rattle(individual, stdev=0.5, x_avg_bond=True):
    """Rattle atoms based on a fraction of the average bond length"""

    if x_avg_bond:
        stdev *= get_avg_radii(individual) * 2

    pos = individual.get_positions()
    individual.set_positions(pos + np.random.normal(scale=stdev, size=pos.shape))

    return
