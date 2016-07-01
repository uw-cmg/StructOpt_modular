import functools
import numpy as np

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
    """A place holder decorator that does nothing. It is purely for documentation."""
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
    if structopt.parameters.globals.USE_MPI4PY:
        from mpi4py import MPI
        # The lists in stuffs_per_core all need to be of the same length 
        amount_of_stuff = 0
        max_stuffs_per_core = max(len(stuffs) for stuffs in stuffs_per_core.values())
        for rank, stuffs in stuffs_per_core.items():
            amount_of_stuff += len(stuffs)
            while len(stuffs) < max_stuffs_per_core:
                stuffs.append(None)

        all_stuffs_per_core = MPI.COMM_WORLD.allgather(stuff)
        correct_stuff = [None for _ in range(amount_of_stuff)]
        for rank, indices in stuffs_per_core.items():
            for index in indices:
                if index is not None:
                    correct_stuff[index] = all_stuffs_per_core[rank][index]

        # If something didn't get sent, use the value on the core
        for i, thing in enumerate(correct_stuff):
            if thing is None:
                correct_stuff[i] = stuff[i]
    else:
        correct_stuff = stuff
    return correct_stuff

