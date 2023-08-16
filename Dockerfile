# FROM ubuntu:latest
# FROM mambaorg/micromamba:jammy
FROM condaforge/mambaforge


ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-base r-cran-randomforest r-cran-devtools python3.6 python3-pip python3-setuptools python3-dev

RUN Rscript -e "devtools::install_github('mjhelf/MassTools')"

WORKDIR GNPS_ML_Processing

COPY /GNPS_ML_Processing/bin/conda_env.yml /conda_env.yml

RUN mamba create -p ./gnps_ml_processing_env
RUN mamba env update -p ./gnps_ml_processing_env --file '/conda_env.yml'