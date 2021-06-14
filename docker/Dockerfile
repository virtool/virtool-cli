FROM continuumio/miniconda3:4.9.2

WORKDIR /app

#add channels
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# install necessary utilities
RUN conda install -c bioconda muscle
RUN conda install -c bioconda cd-hit
RUN conda install -c bioconda hmmer
RUN conda install -c bioconda blast
RUN conda install -c bioconda mcl

COPY . .

RUN pip install .

ENTRYPOINT ["virtool", "vfam"]