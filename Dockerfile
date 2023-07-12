FROM continuumio/miniconda3:23.3.1-0 as base
WORKDIR /app
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda install -c bioconda muscle
RUN conda install -c bioconda cd-hit
RUN conda install -c bioconda hmmer
RUN conda install -c bioconda blast
RUN conda install -c bioconda mcl
RUN pip install poetry
COPY . .

FROM base as build
RUN poetry install --no-dev

FROM base as test
RUN poetry install
RUN poetry run pytest