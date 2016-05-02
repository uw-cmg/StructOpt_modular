def setup(parameter_file):
    import time
    import logging
    import fileio

    parameters = fileio.parameters.read(parameter_file)
    if parameters.globals.USE_MPI4PY:
        try:
            from mpi4py import MPI
        except ImportError:
            raise ImportError("mpi4py must be installed to use StructOpt.")
        parameters.globals.rank = MPI.COMM_WORLD.Get_rank()
    else:
        parameters.globals.rank = 0

    # Setup all the loggers
    logging_level = parameters.globals.get("logging_level", "info")
    logging_level = getattr(logging, logging_level.upper())

    if "loggername" not in parameters.globals:
        parameters.globals.loggername = "{0}-rank{1}-{2}.log".format(parameters.globals.output_filename, parameters.globals.rank, time.strftime("%Y_%m%d_%H%M%S"))
    else:
        raise ValueError("'loggername' should not be defined in the parameter file currently. If you think you want to define it, talk to the developers about why.")

    if parameters.globals.rank == 0:
        logger = fileio.logger_utils.initialize_logger(filename=parameters.globals.output_filename, name="output", level=logging_level)
        default_logger = fileio.logger_utils.initialize_logger(filename=parameters.globals.loggername, name="default", level=logging_level)

        if logging_level <= logging.DEBUG:
            debug_logger = fileio.logger_utils.initialize_logger(filename=parameters.globals.loggername, name="debug", level=logging_level)

    else:
        logger = fileio.logger_utils.initialize_logger(filename=parameters.globals.output_filename, name="output", level=logging_level, disable_output=True)
        default_logger = fileio.logger_utils.initialize_logger(filename=parameters.globals.loggername, name="default", level=logging_level, disable_output=True)

    logger_by_rank = fileio.logger_utils.initialize_logger(filename=parameters.globals.loggername, name="by-rank", level=logging_level)

    if logging_level <= logging.DEBUG:
        debug_logger = fileio.logger_utils.initialize_logger(filename=parameters.globals.loggername, name="by-rank-debug", level=logging_level)

    # Set defaults for parameters and globally scope them
    parameters = fileio.parameters.set_default(parameters)
    globals()["parameters"] = parameters

    structopt.fileio.parameters.write(self)


    return None

