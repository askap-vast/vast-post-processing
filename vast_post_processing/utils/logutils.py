"""Utilities to set up logging. 
"""


# Imports


import logging
import logging.handlers
import datetime, time


# Classes


class UTCFormatter(logging.Formatter):
    """Logging formatter which uses UTC."""

    converter = time.gmtime


# Functions


def setup_logger(verbose: bool, debug: bool, module: str = "full") -> logging.Logger:
    """Set up logging functionality for this module.

    Parameters
    ----------
    verbose : bool
        Flag to display program status and progress to output.
    debug : bool
        Flag to display program errors and actions to output.
    module : str, optional
        Name of output subdirectory indicating which module is being run.
        Defaults to "full", referring to the full core run.

    Returns
    -------
    logging.Logger
        The main Logger object for this module.
    """
    from .. import DATA_SUBDIRECTORIES

    # Log filename as UTC time at run, in ISO-8601
    log_filename = (
        datetime.datetime.now().astimezone().replace(microsecond=0).isoformat()
    )

    # Set up logging level
    logging_level = "WARNING"
    if verbose:
        logging_level = "INFO"
    if debug:
        logging_level = "DEBUG"

    # Create and return logger object
    main_logger = create_logger(
        DATA_SUBDIRECTORIES[module] / log_filename, logging_level
    )
    return main_logger


def create_logger(filename: str, level: str = "WARNING"):
    """
    Create a central file logger with appropriate formatting and dates in UTC.

    Parameters
    ----------
    filename : str
        Filename of output log file.
    level : str, optional
        Logging level to output, by default "WARNING".

    Returns
    -------
    logger
        Logger object
    """
    # Check level is recognized
    level = level.upper()
    assert level in {
        "CRITICAL",
        "ERROR",
        "WARNING",
        "INFO",
        "DEBUG",
    }, "Level must be set to one of the standard Python levels"

    # Add global logging level success
    add_logging_level("SUCCESS", logging.INFO + 5)

    # Add formatter for file
    loggingFormatter = UTCFormatter(
        fmt="%(asctime)-20s %(levelname)-8s %(message)s [%(name)s]",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )

    # Add handler for file
    handler = logging.handlers.FileHandler(filename=filename)
    handler.setFormatter(loggingFormatter)
    handler.setLevel(level)

    # Add formatter for output
    streamloggingFormatter = UTCFormatter(
        fmt=(
            "\033[0;31m%(asctime)-20s\033[39m %(levelname)-8s "
            " \033[0;33m%(message)s\033[39m  \033[0;34m[%(name)s]\033[39m"
        ),
        datefmt="%Y-%m-%dT%H:%M:%S",
    )

    # Add handler for output
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(streamloggingFormatter)
    stream_handler.setLevel(level)

    # Add handlers to logger object and return
    logger = logging.getLogger()
    logger.addHandler(handler)
    logger.addHandler(stream_handler)
    logger.setLevel(level)
    return logger


def add_logging_level(name: str, number: int, method: str = None):
    """Comprehensively add a new logging level to the `logging` module and the
    currently configured logging class.

    Parameters
    ----------
    name : str
        Name of logging level.
    number : int
        Number value of logging level.
    method : str, optional
        Name of convenience method for logging level, defaults to name.lower().

    Notes
    ----------
    Sourced from
    https://stackoverflow.com/a/35804945

    Examples
    --------
    >>> addLoggingLevel('TRACE', logging.DEBUG - 5)
    >>> logging.getLogger(__name__).setLevel("TRACE")
    >>> logging.getLogger(__name__).trace('that worked')
    >>> logging.trace('so did this')
    >>> logging.TRACE
    5

    """
    if not method:
        method = name.lower()

    # This method was inspired by the answers to Stack Overflow post
    # http://stackoverflow.com/q/2183233/2988730, especially
    # http://stackoverflow.com/a/13638084/2988730
    def log_for_level(self, message, *args, **kwargs):
        if self.isEnabledFor(number):
            self._log(number, message, args, **kwargs)

    def log_to_root(message, *args, **kwargs):
        logging.log(number, message, *args, **kwargs)

    logging.addLevelName(number, name)
    setattr(logging, name, number)
    setattr(logging.getLoggerClass(), method, log_for_level)
    setattr(logging, method, log_to_root)
