import time
import logging
import os

from .io.parameters import read as read_parameters
from .io.parameters import write as write_parameters
from .io.logger_utils import initialize_logger, initialize_logger_for_root
from . import common
from .tools import get_size, get_rank

def setup(parameter_file):
    # Read in the parameters
    parameters = read_parameters(parameter_file)

    # Setup all the loggers
    logging_level = logging.parameters.get("logging_level", "info")
    logging_level = getattr(logging, logging_level.upper())

    rank = logging.parameters.rank
    path = logging.parameters.path
    os.makedirs(path, exist_ok=True)

    logger = initialize_logger_for_root(rank=rank, filename=os.path.join(path, 'output.log'), name="output", level=logging_level)
    logger_by_rank = initialize_logger(filename=os.path.join(path, 'log-by-rank-{}.log'.format(rank)), name="by-rank", level=logging_level)

    default_logger = initialize_logger_for_root(rank=rank, filename=os.path.join(path, 'default.log'), name="default", level=logging_level)

    fitness_logger = initialize_logger_for_root(rank=rank, filename=os.path.join(path, 'fitnesses.log'), name="fitness", level=logging_level)

    genealogy_logger = initialize_logger_for_root(rank=rank, filename=os.path.join(path, 'genealogy.log'), name="genealogy", level=logging_level)

    if logging_level <= logging.DEBUG:
        debug_logger = initialize_logger_for_root(rank=rank, filename=os.path.join(path, 'debug.log'), name="debug", level=logging_level)
        debug_logger_by_rank = initialize_logger(filename=os.path.join(path, 'debug-by-rank-{}.log'.format(rank)), name="debug-by-rank", level=logging_level)

    write_parameters(parameters)
    return parameters

__all__ = ['parameters', 'cluster', 'crystal', 'defect', 'surface']
