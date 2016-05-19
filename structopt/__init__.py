import time
import logging
import os

from . import io as structopt_io
from .optimizer import Optimizer
from . import common

def setup(parameter_file):
    parameters = structopt_io.parameters.read(parameter_file)

    # Setup all the loggers
    logging_level = parameters.globals.get("logging_level", "info")
    logging_level = getattr(logging, logging_level.upper())

    if "loggername" not in parameters.globals:
        parameters.globals.loggername = "{0}-rank{1}-{2}".format(parameters.globals.output_filename, parameters.globals.rank, time.strftime("%Y%m%d%H%M%S"))
    else:
        raise ValueError("'loggername' should not be defined in the parameter file currently. If you think you want to define it, talk to the developers about why.")

    if parameters.globals.rank == 0:
        logger = structopt_io.logger_utils.initialize_logger(filename='{}.out'.format(parameters.globals.output_filename), name="output", level=logging_level)
        default_logger = structopt_io.logger_utils.initialize_logger(filename='{}.log'.format(parameters.globals.loggername), name="default", level=logging_level)

        if logging_level <= logging.DEBUG:
            debug_logger = structopt_io.logger_utils.initialize_logger(filename='{}.debug'.format(parameters.globals.loggername), name="debug", level=logging_level)

    else:
        logger = structopt_io.logger_utils.initialize_logger(filename='{}.out'.format(parameters.globals.output_filename), name="output", level=logging_level, disable_output=True)
        default_logger = structopt_io.logger_utils.initialize_logger(filename='{}.log'.format(parameters.globals.loggername), name="default", level=logging_level, disable_output=True)

    logger_by_rank = structopt_io.logger_utils.initialize_logger(filename='{}-by-rank.out'.format(parameters.globals.loggername), name="by-rank", level=logging_level)

    if logging_level <= logging.DEBUG:
        debug_logger = structopt_io.logger_utils.initialize_logger(filename='{}-by-rank.debug'.format(parameters.globals.loggername), name="by-rank-debug", level=logging_level)

    os.makedirs(parameters.globals.output_filename, exist_ok=True)

    # Set defaults for parameters and globally scope them
    parameters = structopt_io.parameters.set_default(parameters)
    globals()["parameters"] = parameters

    structopt_io.parameters.write(parameters)


    return None

__all__ = ['parameters', 'Optimizer', 'cluster', 'crystal', 'defect', 'surface']
