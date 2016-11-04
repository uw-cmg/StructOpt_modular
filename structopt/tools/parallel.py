import sys
import functools
import numpy as np

import structopt


def get_rank():
    if 'mpi4py' in sys.modules:
        from mpi4py import MPI
        return MPI.COMM_WORLD.Get_rank()
    else:
        return 0


def get_size():
    if 'mpi4py' in sys.modules:
        from mpi4py import MPI
        return MPI.COMM_WORLD.Get_size()
    else:
        return 0


def root(method=None, broadcast=True):
    """A decorator to make the function only run on the root node. The returned
    data from the root is then broadcast to all the other nodes and each node
    returns the root's data.
    """
    # This is a complicated decorator. See the following link for some context:
    # https://blogs.it.ox.ac.uk/inapickle/2012/01/05/python-decorators-with-optional-arguments/
    if method is None:
        return functools.partial(root, broadcast=broadcast)

    @functools.wraps(method)
    def wrapper(*args, **kwargs):
        if broadcast and 'mpi4py' in sys.modules:
            from mpi4py import MPI
            if MPI.COMM_WORLD.Get_rank() == 0:
                data = method(*args, **kwargs)
            else:
                data = None
            if hasattr(data, 'bcast'):
                data.bcast()
            else:
                data = MPI.COMM_WORLD.bcast(data, root=0)
        else:
            data = method(*args, **kwargs)
        return data
    wrapper.__doc__ += "\n\n(@root) Designed to run on the root node only.\n"
    return wrapper


def parallel(method):
    """A decorator that does nothing except document that the function is
    designed to run in parallel.
    """
    @functools.wraps(method)
    def wrapper(*args, **kwargs):
        return method(*args, **kwargs)
    wrapper.__doc__ += ("\n\n(@parallel) Designed to run code that runs differently on different cores.\n"
                       "The MPI functionality should be implemented inside these functions.\n")
    return wrapper


def single_core(method):
    """A place holder decorator that does nothing except document that the function is designed to be run on a single core."""
    @functools.wraps(method)
    def wrapper(*args, **kwargs):
        return method(*args, **kwargs)
    return wrapper


def allgather(stuff, stuffs_per_core):
    """Performs an MPI.allgather on a selection of data and uses stuffs_per_core
    to parse out the correct information and return it.

    Args:
        stuff (any): any piece of data (e.g. fitnesses), some of which have been updated
            on their respective cores and some of which haven't. each piece of data should
            be of the same length
        stuffs_per_core (dict<int, list<int>>): a dictionary containing a mapping of the
            cores that contain the correct information to the corresponding indices in the
            pieces of data

    Returns:
        type(stuff): the correct stuff is returned on each core

    Example:
        This is going to take the values::

            values = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

        and convert each of them to strings.

        In this example there are 5 cores, so ``stuffs_per_core`` looks like::

            stuffs_per_core = {0: [0, 5], 1: [1, 6], 2: [2, 7], 3: [3, 8], 4: [4, 9]}

        Now for the code that precedes ``allgather()`` and then calls ``allgather()``::

            # This for-loop modifies different parts of `values` on each core by
            # converting some elements in `values` from an int to a str.
            # We then want to collect the values that each core independently updated
            # and allgather them so that every core has all of the updated values,
            # even though each core only did part of the work.
            for i in stuffs_per_core[rank]:
                values[i] = str(inds[i])
            x = allgather(values, stuffs_per_core)
            print(x)  # returns:  ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
    """
    if hasattr(stuff, 'allgather'):
        raise TypeError('instance of {} has an `allgather` function that should be used instead'.format(stuff.__class__.__name__))

    from mpi4py import MPI
    # The lists in stuffs_per_core all need to be of the same length 
    max_stuffs_per_core = max(len(stuffs) for stuffs in stuffs_per_core.values())
    for rank, stuffs in stuffs_per_core.items():
        while len(stuffs) < max_stuffs_per_core:
            stuffs.append(None)

    all_stuffs_per_core = MPI.COMM_WORLD.allgather(stuff)
    correct_stuff = [None for _ in range(len(stuff))]
    for rank, indices in stuffs_per_core.items():
        for index in indices:
            if index is not None:
                correct_stuff[index] = all_stuffs_per_core[rank][index]

    # If something didn't get sent, use the value on the core
    for i, thing in enumerate(correct_stuff):
        if thing is None:
            correct_stuff[i] = stuff[i]

    return correct_stuff


def parse_MPMD_cores_per_structure(value):
    """Converts an input ``value`` from a value in the parameter file into a ``{'min': ..., 'max': ...}`` dictionary."""
    if isinstance(value, int):
        if int == 0:
            return None
        else:
            return {'min': value, 'max': value}
    elif isinstance(value, str):
        if value == 'any':
            from mpi4py import MPI
            return {'min': 1, 'max': MPI.COMM_WORLD.Get_size()}
        elif '-' not in value:
            value = int(value)
            return {'min': value, 'max': value}
        else:
            min_, max_ = value.split('-')
            return {'min': int(min_), 'max': int(max_)}
    else:
        raise TyeError("'MPMD_cores_per_structure' must be an 'int' or 'str'")


class SingleCorePerIndividual(object):
    """
    This is an untested context manager for using the mpi4py parallelization technique.
    I'm not sure we should use it, but I'm thinking on it.
    The nice thing is that it is self-contained reasonably well.
    The bad thing is that the user has to manually set fields on the individuals before
    the context manager exits in order for it to work properly,
    and there isn't a great way to error check whether they did that or not.

    Example usage:

        def calculate_energies(population):
            to_fit = [individual for individual in population if individual.LAMMPS is None]
            with SingleCorePerIndividual(population, to_run, "LAMMPS") as to_compute:
                for individual in to_compute:
                    individual.LAMMPS = individual.fitnesses.LAMMPS.fitness(individual)
            return [individual.LAMMPS for individual in population]
    """
    def __init__(population, to_run, attr_name):
        self.population = population
        self.to_run = to_run
        import parameters
        self.rank = parameters.rank
        self.ncores = parameters.ncores

    def __enter__(self):
        self.individuals_per_core = {r: [] for r in range(self.ncores)}
        for i, individual in enumerate(self.to_run):
            self.individuals_per_core[i % self.ncores].append(individual)
        return individuals_per_core[rank]

    def __exit__(self, exception_type, exception_value, traceback):
        positions_per_core = {rank: [self.population.position(individual) for individual in individuals] for rank, individuals in self.individuals_per_core.items()}
        self.allgather([getattr(individual, self.attr_name) for individual in self.population])
        allgather(results, positions_per_core)

