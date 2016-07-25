import time
import logging
import os

from . import io as structopt_io
from . import common
from .tools import get_size, get_rank

def setup(parameter_file):
    parameters = structopt_io.parameters.read(parameter_file)

    # Setup all the loggers
    logging_level = logging.parameters.get("logging_level", "info")
    logging_level = getattr(logging, logging_level.upper())

    if "loggername" not in logging.parameters:
        logging.parameters.loggername = "{0}-rank{1}-{2}".format(logging.parameters.output_filename, logging.parameters.rank, time.strftime("%Y%m%d%H%M%S"))
    else:
        raise ValueError("'loggername' should not be defined in the parameter file currently. If you think you want to define it, talk to the developers about why.")

    if logging.parameters.rank == 0:
        logger = structopt_io.logger_utils.initialize_logger(filename='{}.out'.format(logging.parameters.output_filename), name="output", level=logging_level)
        default_logger = structopt_io.logger_utils.initialize_logger(filename='{}.log'.format(logging.parameters.loggername), name="default", level=logging_level)

        if logging_level <= logging.DEBUG:
            debug_logger = structopt_io.logger_utils.initialize_logger(filename='{}.debug'.format(logging.parameters.loggername), name="debug", level=logging_level)

    else:
        logger = structopt_io.logger_utils.initialize_logger(filename='{}.out'.format(logging.parameters.output_filename), name="output", level=logging_level, disable_output=True)
        default_logger = structopt_io.logger_utils.initialize_logger(filename='{}.log'.format(logging.parameters.loggername), name="default", level=logging_level, disable_output=True)

    logger_by_rank = structopt_io.logger_utils.initialize_logger(filename='{}-by-rank.out'.format(logging.parameters.loggername), name="by-rank", level=logging_level)

    if logging_level <= logging.DEBUG:
        debug_logger = structopt_io.logger_utils.initialize_logger(filename='{}-by-rank.debug'.format(logging.parameters.loggername), name="by-rank-debug", level=logging_level)

    os.makedirs(logging.parameters.output_filename, exist_ok=True)

    structopt_io.parameters.write(parameters)

    return parameters

__all__ = ['parameters', 'cluster', 'crystal', 'defect', 'surface']
