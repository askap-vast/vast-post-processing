"""Utilities to set up logging. 
"""


# Imports


import logging
import logging.handlers
import time


# Constants


# Setting Logger Rotation Settings
LOGGER_MAX_BYTES = 1000000
"""Maximum number of bytes for the log to contain.
Rollover occurs when this amount is reached in one log.
"""


LOGGER_BACKUP_COUNT = 4
"""Maximum number of backup log files. 
When rollover occurs, new files with number suffices up to this amount will be
created.
"""


# Classes


class UTCFormatter(logging.Formatter):
    """Logging formatter which uses UTC."""

    converter = time.gmtime


# Functions


def create_logger(filename: str, level: str = "WARNING"):
    """
    Create a central rotating file logger with the appropriate
    formatting and dates in UTC

    Parameters
    ----------
    filename : str
        Filename of output log file
    level : str, optional
        Logging level to output, by default "WARNING"

    Returns
    -------
    logger
        Logger object
    """

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

    loggingFormatter = UTCFormatter(
        fmt="%(asctime)-20s %(levelname)-8s %(message)s [%(name)s]",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )

    handler = logging.handlers.RotatingFileHandler(
        filename=filename, maxBytes=LOGGER_MAX_BYTES, backupCount=LOGGER_BACKUP_COUNT
    )
    handler.setFormatter(loggingFormatter)
    handler.setLevel(level)

    # Adding Stream Handler

    streamloggingFormatter = UTCFormatter(
        fmt=(
            "\033[0;31m%(asctime)-20s\033[39m %(levelname)-8s "
            " \033[0;33m%(message)s\033[39m  \033[0;34m[%(name)s]\033[39m"
        ),
        datefmt="%Y-%m-%dT%H:%M:%S",
    )

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(streamloggingFormatter)
    stream_handler.setLevel(level)

    # Adding to Logger

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
