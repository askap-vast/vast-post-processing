ARG VIRTUAL_ENV=/venv

FROM pawsey/hpc-python:2021.09 as base
ENV PYTHONUNBUFFERED=1
WORKDIR /app

FROM base as builder
ARG VIRTUAL_ENV
ENV POETRY_VERSION=1.1.8 \
    PIP_NO_CACHE_DIR=1 \
    VIRTUAL_ENV=$VIRTUAL_ENV
RUN pip install "poetry==${POETRY_VERSION}"
ENV DEBIAN_FRONTEND=noninteractive
COPY pyproject.toml poetry.lock ./
RUN poetry export -f requirements.txt | pip install -r /dev/stdin

FROM builder as final
COPY link_neighbours.py correct_vast.py convolve_neighbours.py ./
