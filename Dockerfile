ARG VIRTUAL_ENV=/venv

FROM pawsey/hpc-python:2021.09 as base_pawsey
ARG VIRTUAL_ENV
ENV PYTHONUNBUFFERED=1 \
    DEBIAN_FRONTEND=noninteractive \
    POETRY_VERSION=1.1.8 \
    PIP_NO_CACHE_DIR=1 \
    VIRTUAL_ENV=$VIRTUAL_ENV
WORKDIR /app
RUN pip install "poetry==${POETRY_VERSION}" && \
    apt-get update && \
    apt-get -y install git swarp && \
    rm -rf /var/lib/apt/lists/*
COPY pyproject.toml poetry.lock ./

FROM python:3.9-bullseye as base_mortimer
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update -qq && \
    apt-get -y --no-install-recommends install \
        build-essential \
        ca-certificates \
        gdb \
        gfortran \
        wget && \
    apt-get clean all && \
    rm -r /var/lib/apt/lists/*
ARG VIRTUAL_ENV=/venv
ARG MPICH_DIR="/opt/mpich"
ARG MPICH_VERSION="3.2.1"

ENV FFLAGS="-w -fallow-argument-mismatch -O2"
RUN mkdir -p /tmp/mpich-build && \
    cd /tmp/mpich-build && \
    wget -O mpich-$MPICH_VERSION.tar.gz http://www.mpich.org/static/downloads/$MPICH_VERSION/mpich-$MPICH_VERSION.tar.gz && \
    tar xzf mpich-$MPICH_VERSION.tar.gz && \
    cd /tmp/mpich-build/mpich-$MPICH_VERSION && \
    ./configure --prefix=$MPICH_DIR && \
    make install && \
    # ldconfig && \
    cd / && \
    rm -rf /tmp/mpich-build

ENV PYTHONUNBUFFERED=1 \
    POETRY_VERSION=1.1.8 \
    POETRY_HOME=/etc/poetry \
    PIP_NO_CACHE_DIR=1 \
    VIRTUAL_ENV=$VIRTUAL_ENV \
    PATH=${MPICH_DIR}/bin:$PATH \
    LD_LIBRARY_PATH=${MPICH_DIR}/lib:$LD_LIBRARY_PATH \
    MANPATH=${MPICH_DIR}/share/man:$MANPATH
RUN mkdir -p ${VIRTUAL_ENV} ${POETRY_HOME} && \
    apt-get update -qq && \
    # apt-get -y install software-properties-common && \
    # add-apt-repository ppa:deadsnakes/ppa && \
    # apt-get update -qq && \
    # apt-get -y install python3.9 python3.9-venv curl git swarp && \
    apt-get -y install curl git swarp && \
    curl -sSL https://install.python-poetry.org | python3 - && \
    apt-get clean all && \
    rm -rf /var/lib/apt/lists/*
ENV PATH=$PATH:${POETRY_HOME}/bin
WORKDIR /app
COPY pyproject.toml poetry.lock ./

FROM base_mortimer as dev
# install python deps, include dev deps
#RUN poetry export -f requirements.txt -E mpi --dev --without-hashes | pip install -r /dev/stdin
RUN poetry config virtualenvs.create false
RUN poetry install -E mpi
COPY link_neighbours.py correct_vast.py convolve_neighbours.py mpi_logger.py swarp.py selavy_combined.py cleanup.py mpi_test.py ./

FROM base_pawsey as final
# install python deps
RUN poetry export -f requirements.txt --without-hashes | pip install -r /dev/stdin
COPY link_neighbours.py correct_vast.py convolve_neighbours.py mpi_logger.py swarp.py selavy_combined.py cleanup.py mpi_test.py ./

FROM base_mortimer as mortimer
RUN poetry export -f requirements.txt -E mpi --without-hashes | pip install -r /dev/stdin
COPY link_neighbours.py correct_vast.py convolve_neighbours.py mpi_logger.py swarp.py selavy_combined.py cleanup.py mpi_test.py ./
