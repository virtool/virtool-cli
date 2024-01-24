FROM continuumio/miniconda3:latest as base
WORKDIR /app
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda install -c bioconda cd-hit hmmer blast mcl muscle -q -y
RUN pip install poetry
COPY pyproject.toml poetry.lock ./
RUN poetry install --no-root

FROM base as test
COPY virtool_cli ./virtool_cli
COPY assets ./assets
COPY tests ./tests
RUN poetry install
ENTRYPOINT ["poetry", "run"]