import logging


def initialize_logger(filename="default.log", name="default", level=logging.INFO, disable_output=False):
    logger = logging.getLogger(name)
    handler = logging.FileHandler(filename)
    formatter = logging.Formatter("%(asctime)s : %(levelname)s : %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(level)

    if disable_output:
        return_none = lambda *args, **kwargs: None
        logger.debug = return_none
        logger.info = return_none
        logger.warning = return_none
        logger.error = return_none
        logger.critical = return_none

    return logger
