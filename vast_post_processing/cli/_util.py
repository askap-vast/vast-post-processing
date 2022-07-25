import multiprocessing
import sys
from typing import Union

from loguru import logger  # noqa: F401
from schwimmbad import choose_pool, SerialPool, MultiPool, MPIPool


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
    global logger
    logger.remove()
    logger.add(
        sys.stderr,
        level=level,
        enqueue=True,
        format=(
            "<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | <level>{level: <8}</level> |"
            " {extra[worker_name]: <18} |"
            " <cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan>"
            " - <level>{message}</level>"
        ),
    )
    logger.configure(extra={"worker_name": _get_worker_name(mpi=mpi, n_proc=n_proc)})

    pool = choose_pool(mpi=mpi, processes=n_proc)
    return pool
