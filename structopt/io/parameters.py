"""Contains functionality for reading, writing, and parsing StrcutOpt parameters."""

import json
import logging
import os
from structopt.tools.dictionaryobject import DictionaryObject


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
    for module in parameters.relaxations.modules:
        parameters.relaxations[module].setdefault('use_mpi4py', False)
        parameters.relaxations[module].setdefault('MPMD_cores_per_structure', 0)
        if parameters.relaxations[module].use_mpi4py:
            use_mpi4py = True
    for module in parameters.fitnesses.modules:
        parameters.fitnesses[module].setdefault('use_mpi4py', False)
        parameters.fitnesses[module].setdefault('MPMD_cores_per_structure', 0)
        if parameters.fitnesses[module].use_mpi4py:
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
    logger = logging.getLogger('default')

    parameters.setdefault('seed', None)

    parameters.logging.setdefault('output_filename', 'Output')

    if 'relaxations' not in parameters or not parameters['relaxations']:
        raise ValueError('Relaxations must be specified in the parameter file.')

    if 'fitnesses' not in parameters or not parameters['fitnesses']:
        raise ValueError('Fitnesses must be specified in the parameter file.')

    parameters.setdefault('weights', [1.0 for _ in parameters.fitnesses.weights])

    parameters.selections.crossover_probability = sum(parameters.crossovers.options.values())

    parameters.convergence.setdefault('maxgen', 10)

    return parameters

