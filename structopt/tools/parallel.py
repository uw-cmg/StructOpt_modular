import functools

import structopt


def root(func):
    """A decorator to make the function only run on the root node. The returned
    data from the root is then broadcast to all the other nodes and each node
    returns the root's data.
    """
    if structopt.parameters.globals.USE_MPI4PY:
        from mpi4py import MPI

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        if structopt.parameters.globals.rank == 0:
            data = func(*args, **kwargs)
        else:
            data = None
        if structopt.parameters.globals.USE_MPI4PY:
            data = MPI.COMM_WORLD.bcast(data, root=0)
        return data
    wrapper.__doc__ += "\nDesigned to run on the root node only.\n"
    return wrapper


def parallel(func):
    """A decorator that does nothing except document that the function is
    designed to run in parallel.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs):
    wrapper.__doc__ += "\nDesigned to run code that runs differently on different cores.\n"
                       "The MPI functionality should be implemented inside these functions.\n"
    return wrapper


def single_core(func):
    """A place holder decorator that does nothing. It is purely for documentation."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs):
    return wrapper

