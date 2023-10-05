"""Utilities for the CLI subpackage.
"""


# Imports


import sys
import logging
import multiprocessing
from typing import Union

from schwimmbad import choose_pool, SerialPool, MultiPool, MPIPool


# Constants


logger = logging.getLogger(__name__)
"""Global reference to the logger for this project.
"""


# Functions


def _get_worker_name(mpi: bool = False, n_proc: int = 1) -> str:
    if mpi:
        from mpi4py import MPI

        worker_name = f"MPI-{MPI.COMM_WORLD.Get_rank()}"
    elif n_proc > 1:
        worker_name = multiprocessing.current_process().name
    else:
        worker_name = "SerialWorker"
    return worker_name


def get_pool(
    level: str = "DEBUG", mpi: bool = False, n_proc: int = 1
) -> Union[SerialPool, MultiPool, MPIPool]:
    # Temporarily reset logger format
    previous_fmt = logger.handlers[0].formatter._fmt
    logger.handlers[0].formatter._fmt = (
        "<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | <level>{level: <8}</level> |"
        " {extra[worker_name]: <18} |"
        " <cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan>"
        " - <level>{message}</level>"
    )

    # Add worker name as extra
    logger = logging.LoggerAdapter(
        extra={"worker_name": _get_worker_name(mpi=mpi, n_proc=n_proc)}
    )

    # Choose pool
    pool = choose_pool(mpi=mpi, processes=n_proc)

    # Reset logger format
    logger.handlers[0].formatter._fmt = previous_fmt
    return pool
