"""Contains functionality for reading, writing, and parsing StrcutOpt parameters."""

import json
import logging
import os
from structopt.tools.dictionaryobject import DictionaryObject
import time
import distutils.spawn

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

    # Make sure we have a logging dictionary
    parameters.setdefault('logging', DictionaryObject({}))

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
            import mpi4py
        except ImportError:
            raise ImportError("mpi4py must be installed to use StructOpt.")
        mpiexec_path, _ = os.path.split(distutils.spawn.find_executable("mpiexec"))
        for executable, path in mpi4py.get_config().items():
            if executable not in ['mpicc', 'mpicxx', 'mpif77', 'mpif90', 'mpifort']:
                continue
            if mpiexec_path not in path:
                raise ImportError("mpi4py may not be configured against the same version of 'mpiexec' that you are using. The 'mpiexec' path is {mpiexec_path} and mpi4py.get_config() returns:\n{mpi4py_config}\n".format(mpiexec_path=mpiexec_path, mpi4py_config=mpi4py.get_config()))
        from mpi4py import MPI
        if 'Open MPI' not in MPI.get_vendor():
            raise ImportError("mpi4py must have been installed against Open MPI in order for StructOpt to function correctly.")
        vendor_number = ".".join([str(x) for x in MPI.get_vendor()[1]])
        if vendor_number not in mpiexec_path:
            raise ImportError("The MPI version that mpi4py was compiled against does not match the version of 'mpiexec'. mpi4py's version number is {}, and mpiexec's path is {}".format(MPI.get_vendor(), mpiexec_path))
        #print("Current MPI configuration:\n")
        #print("'mpiexec' path:")
        #print(mpiexec_path)
        #print("mpi4py.get_config():")
        #print(mpi4py.get_config())
        #print("MPI.get_vendor():")
        #print(MPI.get_vendor())


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
    logging.parameters.generation = 0

    if "path" in logging.parameters:
        raise ValueError("'path' should not be defined in the parameter file currently. If you think you want to define it, talk to the developers about why.")
    
    # If parallel and no seed, all nodes need the same seed and same path
    if parameters.logging.ncores > 1:
        from mpi4py import MPI
        seed = MPI.COMM_WORLD.bcast(int(time.time()), root=0)
        path = MPI.COMM_WORLD.bcast(os.path.join(os.getcwd(), "logs{}".format(time.strftime("%Y%m%d%H%M%S"))), root=0)
    else:
        seed = None
        path = os.path.join(os.getcwd(), "logs{}".format(time.strftime("%Y%m%d%H%M%S")))

    logging.parameters.path = path        
    parameters.setdefault('seed', seed)
    parameters.setdefault('post_processing', DictionaryObject({}))
    parameters.post_processing.setdefault('XYZs', 0)
    parameters.setdefault('adaptation', [])
    parameters.setdefault('fingerprinters', DictionaryObject({'options': []}))
    parameters.convergence.setdefault('max_generations', 10)

    if 'relaxations' not in parameters or not parameters['relaxations']:
        raise ValueError('Relaxations must be specified in the parameter file.')

    if 'fitnesses' not in parameters or not parameters['fitnesses']:
        raise ValueError('Fitnesses must be specified in the parameter file.')

    # Make sure every operation is defined, and that every operation has a
    # kwargs
    for operation in MODULES:
        parameters.setdefault(operation, None)
        if parameters[operation] is not None:
            for operator in parameters[operation]:
                parameters[operation][operator].setdefault('kwargs', {})

    return parameters

