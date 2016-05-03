import functools

import structopt
from structopt.common.individual import Individual


def generate(individual, *args, **kwargs):
    """ Uses the relevant parameters from structopt to intialize the input Individual by modifying it in-place.

        Args:
            individual (Individual): an Individual that is uninitialized
            *args, **kwargs: arguments for either ase.Atoms or a different generator function
    """
    return None

