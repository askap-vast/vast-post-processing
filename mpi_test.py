import logging
import os
import socket
from typing import Any

from mpi4py import MPI
import schwimmbad
import structlog
import typer

import vast_post_processing.mpi_logger as mpi_logger

slurm_job_id = os.environ.get("SLURM_JOB_ID", "no-slurm")

# configure root logger to use structlog
structlog.configure(
    processors=mpi_logger.LOGGING_COMMON_PROCESSORS,  # type: ignore
    logger_factory=structlog.stdlib.LoggerFactory(),
)
HANDLER = mpi_logger.MPIFileHandler(f"swarp-{slurm_job_id}.log")
FORMATTER = logging.Formatter("%(message)s")
HANDLER.setFormatter(FORMATTER)

LOGGER = logging.getLogger("mpi_test")
LOGGER.setLevel(logging.DEBUG)
LOGGER.addHandler(HANDLER)


def test_worker(args: Any):
    _logger = structlog.get_logger("mpi_test")
    logger = _logger.bind(rank=MPI.COMM_WORLD.Get_rank(), hostname=socket.gethostname())
    logger.debug(f"worker args: {args}")


def main(n_inputs: int = 103, n_proc: int = 1, mpi: bool = False):
    logger = structlog.get_logger("mpi_test").bind(
        rank=MPI.COMM_WORLD.Get_rank(),
        size=MPI.COMM_WORLD.Get_size(),
        hostname=socket.gethostname(),
    )
    logger.info("checking rank and size")
    pool = schwimmbad.choose_pool(mpi=mpi, processes=n_proc)

    # if using MPI, the following is executed only on the main process
    arg_list: list[str] = []
    for i in range(n_inputs):
        arg_list.append(f"Input {i:03d}")
        logger.info(f"Added input to arg list: {i:03d}")

    # distribute tasks
    pool.map(test_worker, arg_list)
    pool.close()


if __name__ == "__main__":
    typer.run(main)
