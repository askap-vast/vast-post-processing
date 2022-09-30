ARG VIRTUAL_ENV="/venv"

FROM python:3.9-bullseye as base_mortimer
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update -qq && \
    apt-get -y --no-install-recommends install \
        build-essential \
        ca-certificates \
        gdb \
        gfortran \
        wget \
        curl \
        git \
        swarp && \
    apt-get clean all && \
    rm -r /var/lib/apt/lists/*
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
WORKDIR /app

FROM base_mortimer as builder_mortimer
ARG VIRTUAL_ENV
ENV PYTHONUNBUFFERED=1 \
    POETRY_VERSION=1.1.8 \
    POETRY_HOME=/etc/poetry \
    # POETRY_VIRTUALENVS_PATH=/venv \
    VIRTUAL_ENV=${VIRTUAL_ENV} \
    PIP_NO_CACHE_DIR=1 \
    PATH=${MPICH_DIR}/bin:$PATH \
    LD_LIBRARY_PATH=${MPICH_DIR}/lib:$LD_LIBRARY_PATH \
    MANPATH=${MPICH_DIR}/share/man:$MANPATH

RUN mkdir -p ${VIRTUAL_ENV} ${POETRY_HOME} && \
    curl -sSL https://install.python-poetry.org | python3 - && \
    python3 -m venv ${VIRTUAL_ENV}
ENV PATH="${VIRTUAL_ENV}/bin:$PATH:${POETRY_HOME}/bin"
RUN poetry config virtualenvs.create false
COPY pyproject.toml poetry.lock ./
RUN poetry install --no-root -E mpi
COPY . .

FROM base_mortimer as dev
ARG VIRTUAL_ENV
ENV VIRTUAL_ENV=${VIRTUAL_ENV} \
    PATH="${VIRTUAL_ENV}/bin:$PATH"
COPY --from=builder_mortimer ${VIRTUAL_ENV} ${VIRTUAL_ENV}

COPY . .
RUN pip install --no-deps .
