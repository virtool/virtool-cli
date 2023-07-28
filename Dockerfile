FROM continuumio/miniconda3:latest as base
WORKDIR /app
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda install -c bioconda cd-hit hmmer blast mcl muscle -q -y
RUN pip install poetry
COPY . .

FROM base as build
RUN poetry install --only main

FROM base as test
RUN poetry install
RUN poetry run pytest