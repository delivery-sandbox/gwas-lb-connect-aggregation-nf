FROM continuumio/miniconda3:4.10.3

################## METADATA ######################
LABEL name="quay.io/lifebitaiorg/harmonisation-mungesumstats" \
      description="Docker image with MungeSumstats and associated R packages" \
      version="1.1.0" \
      maintainer="Ilya Levantis <ilya@lifebit.ai>"

################## INSTALLATION ######################
USER root
ARG ENV_NAME="base"

RUN apt-get --allow-releaseinfo-change update -y \
    && apt-get install -y \
              procps \
              wget \
              libcurl4-openssl-dev \
              build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install mamba for faster installation in the subsequent step
RUN conda install -c conda-forge mamba -y

# Install the conda environment
COPY environment.yml /
RUN mamba env update --quiet --name ${ENV_NAME} --file /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/${ENV_NAME}/bin:$PATH

COPY install.R /
RUN Rscript /install.R

# Dump the details of the installed packages to a file for posterity
RUN mamba env export --name ${ENV_NAME} > ${ENV_NAME}_exported.yml

CMD ["bash"]
