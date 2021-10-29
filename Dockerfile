ARG VIRTUAL_ENV=/venv

FROM pawsey/hpc-python:2021.09 as base
ENV PYTHONUNBUFFERED=1
WORKDIR /app

FROM base as builder
ARG VIRTUAL_ENV
ENV POETRY_VERSION=1.1.8 \
    PIP_NO_CACHE_DIR=1 \
    VIRTUAL_ENV=$VIRTUAL_ENV
RUN pip install "poetry==${POETRY_VERSION}" && \
    apt-get update && \
    apt-get -y install git swarp && \
    rm -rf /var/lib/apt/lists/*
ENV DEBIAN_FRONTEND=noninteractive
COPY pyproject.toml poetry.lock ./

FROM builder as dev
# install python deps, include dev deps
RUN poetry export -f requirements.txt --dev --without-hashes | pip install -r /dev/stdin
COPY link_neighbours.py correct_vast.py convolve_neighbours.py swarp.py selavy_combined.py cleanup.py ./

FROM builder as final
# install python deps
RUN poetry export -f requirements.txt --without-hashes | pip install -r /dev/stdin
COPY link_neighbours.py correct_vast.py convolve_neighbours.py swarp.py selavy_combined.py cleanup.py ./
