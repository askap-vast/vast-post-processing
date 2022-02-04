FROM pawsey/hpc-python:2021.09 as base
ENV PYTHONFAULTHANDLER=1 \
    PYTHONHASHSEED=random \
    PYTHONUNBUFFERED=1
WORKDIR /app

FROM base as builder-base
ENV POETRY_VERSION=1.1.12 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    DEBIAN_FRONTEND=noninteractive
RUN pip install "poetry==${POETRY_VERSION}" && \
    python -m venv /venv && \
    apt-get update && \
    apt-get -y install git swarp && \
    rm -rf /var/lib/apt/lists/*
COPY poetry.lock pyproject.toml ./

FROM builder-base as builder-prod
RUN poetry export -f requirements.txt --without-hashes --no-interaction --no-ansi | /venv/bin/pip install -r /dev/stdin
COPY . .
RUN poetry build && /venv/bin/pip install dist/*.whl

FROM builder-base as builder-dev
RUN poetry export -f requirements.txt --without-hashes --dev --no-interaction --no-ansi | /venv/bin/pip install -r /dev/stdin
COPY . .
RUN poetry build && /venv/bin/pip install dist/*.whl

FROM base as prod
COPY --from=builder-prod /venv /venv
ENV PATH="/venv/bin:${PATH}" \
    VIRTUAL_ENV="/venv"

FROM builder-dev as dev
ENV PATH="/venv/bin:${PATH}" \
    VIRTUAL_ENV="/venv"
