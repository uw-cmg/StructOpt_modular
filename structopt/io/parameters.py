"""Contains functionality for reading, writing, and parsing StrcutOpt parameters."""

import json
import logging
import os
from structopt.tools.dictionaryobject import DictionaryObject
import time

MODULES = ['relaxations', 'fitnesses', 'mutations', 'generators', 'crossovers', 'selections', 'predators', 'pso_moves']


def read(input):
    """Sets StructOpt parameters from a dictionary or filename"""

    # Read and store the input
    if isinstance(input, dict):
        parameters = input
    elif isinstance(input, str) and os.path.exists(input):
        parameters = DictionaryObject(json.load(open(input)))
    else:
        raise IOError('Error in input or input file')

    # If mpi4py is used, make sure we can import it and set the rank/size for all cores in the logging parameters
    use_mpi4py = False
    for module in parameters.relaxations:
        parameters.relaxations[module].kwargs.setdefault('use_mpi4py', False)
        parameters.relaxations[module].kwargs.setdefault('MPMD_cores_per_structure', 0)
        if parameters.relaxations[module].kwargs.use_mpi4py:
            use_mpi4py = True
    for module in parameters.fitnesses:
        parameters.fitnesses[module].kwargs.setdefault('use_mpi4py', False)
        parameters.fitnesses[module].kwargs.setdefault('MPMD_cores_per_structure', 0)
        if parameters.fitnesses[module].kwargs.use_mpi4py:
            use_mpi4py = True

    if use_mpi4py:
        try:
            from mpi4py import MPI
        except ImportError:
            raise ImportError("mpi4py must be installed to use StructOpt.")
        if 'Open MPI' not in MPI.get_vendor():
            raise ImportError("mpi4py must have been installed against Open MPI in order for StructOpt to function correctly.")
        parameters.logging.rank = MPI.COMM_WORLD.Get_rank()
        parameters.logging.ncores = MPI.COMM_WORLD.Get_size()
    else:
        parameters.logging.rank = 0
        parameters.logging.ncores = 1

    # Save the logging parameters in the actual logging module
    logging.parameters = parameters.logging

    parameters = set_default(parameters)
    return parameters


def write(parameters):
    output = logging.getLogger('output')
    output.info('Current parameters:')
    output.info(json.dumps(parameters, sort_keys=True, indent=4))


def set_default(parameters):
    if "path" not in logging.parameters:
        logging.parameters.path = os.path.join(os.getcwd(), "logs{}".format(time.strftime("%Y%m%d%H%M%S")))
    else:
        raise ValueError("'path' should not be defined in the parameter file currently. If you think you want to define it, talk to the developers about why.")
    logging.parameters.generation = 0

    # If parallel and no seed, all nodes need the same seed
    if parameters.logging.ncores > 1:
        from mpi4py import MPI
        seed = MPI.COMM_WORLD.bcast(int(time.time()), root=0)
    else:
        seed = None

    parameters.setdefault('seed', seed)
    parameters.setdefault('post_processing', DictionaryObject({}))
    parameters.post_processing.setdefault('XYZs', 0)
    parameters.setdefault('adaptation', [])

    if 'relaxations' not in parameters or not parameters['relaxations']:
        raise ValueError('Relaxations must be specified in the parameter file.')

    if 'fitnesses' not in parameters or not parameters['fitnesses']:
        raise ValueError('Fitnesses must be specified in the parameter file.')

    parameters.convergence.setdefault('max_generations', 10)
    for module_name in MODULES:
        parameters.setdefault(module_name, None)

    # Make sure every operation has a kwargs. Not sure about fingerprinters yet.
    for operation in MODULES:
        if parameters[operation] is None:
            continue
        for operator in parameters[operation]:
            parameters[operation][operator].setdefault('kwargs', {})

    return parameters

