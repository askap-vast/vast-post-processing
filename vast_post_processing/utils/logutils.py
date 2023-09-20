"""
Logging Setup Functions
"""

import logging
import logging.handlers
import time

# Setting Logger Rotation Settings
LOGGER_MAX_BYTES = 1000000
LOGGER_BACKUP_COUNT = 4


class UTCFormatter(logging.Formatter):
    """
    Creating a logging formatter that uses UTC
    """

    converter = time.gmtime


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