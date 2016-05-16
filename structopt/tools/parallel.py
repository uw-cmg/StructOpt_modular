import functools

import structopt


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
        if structopt.parameters.globals.rank == 0:
            data = method(*args, **kwargs)
        else:
            data = None
        if structopt.parameters.globals.USE_MPI4PY:  # This if statement exists because the code shouldn't break when not using mpi4py
            from mpi4py import MPI
            data = MPI.COMM_WORLD.bcast(data, root=0)
        return data
    wrapper.__doc__ += "\nDesigned to run on the root node only.\n"
    return wrapper


def parallel(method):
    """A decorator that does nothing except document that the function is
    designed to run in parallel.
    """
    @functools.wraps(method)
    def wrapper(*args, **kwargs):
        return method(*args, **kwargs)
    wrapper.__doc__ += ("\nDesigned to run code that runs differently on different cores.\n"
                       "The MPI functionality should be implemented inside these functions.\n")
    return wrapper


def single_core(method):
    """A place holder decorator that does nothing. It is purely for documentation."""
    @functools.wraps(method)
    def wrapper(*args, **kwargs):
        return method(*args, **kwargs)
    return wrapper

