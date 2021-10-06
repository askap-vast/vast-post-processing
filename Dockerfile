ARG VIRTUAL_ENV=/venv

FROM python:3.9-buster as base
ENV PYTHONUNBUFFERED=1
WORKDIR /app

FROM base as builder
ARG VIRTUAL_ENV
ENV POETRY_VERSION=1.1.8 \
    PIP_NO_CACHE_DIR=1 \
    VIRTUAL_ENV=$VIRTUAL_ENV
RUN pip install "poetry==${POETRY_VERSION}"
RUN python -m venv ${VIRTUAL_ENV}
ENV PATH="${VIRTUAL_ENV}/bin:${PATH}"
ENV DEBIAN_FRONTEND=noninteractive

# RACS-tools installation requires numpy and a Fortran compiler
RUN apt-get update && \
    apt-get -y install gfortran && \
    rm -rf /var/lib/apt/lists/* && \
    pip install numpy
COPY pyproject.toml poetry.lock ./
RUN poetry export -f requirements.txt | pip install -r /dev/stdin

FROM base as final
ARG VIRTUAL_ENV
ENV PATH="${VIRTUAL_ENV}/bin:${PATH}"
RUN apt-get update && \
    apt-get -y install libgfortran5 && \
    rm -rf /var/lib/apt/lists/*
COPY --from=builder $VIRTUAL_ENV $VIRTUAL_ENV
COPY link_neighbours.py correct_vast.py convolve_neighbours.py ./
