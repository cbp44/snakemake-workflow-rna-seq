############ miniconda3 ####################
FROM docker.io/continuumio/miniconda3:4.10.3 AS miniconda3

COPY ./.docker/condarc /opt/conda/.condarc

RUN set -eux; \
  # Install mamba for faster dependency resolution than conda
  conda install -y -c conda-forge \
    conda-forge::mamba==0.19.1; \
  # Clean up anything left over to trim down the image
  conda clean --all -y; \
  # Replace http with https in apt config
  sed -i 's/^deb http:/deb https:/g' /etc/apt/sources.list; \
  # Add shell init code to activate conda and load its base environment
  conda init bash; \
  echo "conda activate base" >> ~/.bashrc; \
  # Clean some files out
	rm -rf /var/lib/apt/lists/* /usr/share/doc/* /usr/share/man/* /opt/conda/man/*
#################### /miniconda3 ####################


#################### snakefmt ####################
FROM miniconda3 AS snakefmt

# Install snakefmt into snakefmt conda environment
RUN set -eu; \
  . /root/.bashrc; \
  mamba create -y \
    -n snakefmt; \
  conda activate snakefmt; \
  pip install --no-cache-dir snakefmt
#################### /snakefmt ####################


#################### snakemake ####################
FROM snakefmt AS snakemake

# Where snakemake conda environments get created
ENV SNAKEMAKE_CONDA_PREFIX=/mnt/snakemake-envs

# Location of snakemake output cache
ENV SNAKEMAKE_OUTPUT_CACHE=/mnt/snakemake-cache

# Define SNAKEMAKE_VERSION build arg that persists in Docker environment
ARG SNAKEMAKE_VERSION
ENV SNAKEMAKE_VERSION=${SNAKEMAKE_VERSION:-6.12.3}

RUN set -eux; \
  mamba create -y \
    -c conda-forge -c bioconda \
    -n snakemake \
    bioconda::snakemake-minimal=="${SNAKEMAKE_VERSION}" \
    conda-forge::singularity==3.8.5; \    
  mamba clean --all -y; \
  rm -rf /opt/conda/man/*; \
  sed -i 's/conda activate base/conda activate snakemake/' ~/.bashrc; \
  mkdir -p \
    "$SNAKEMAKE_CONDA_PREFIX" "SNAKEMAKE_OUTPUT_CACHE"; \
  echo "export SNAKEMAKE_OUTPUT_CACHE=$SNAKEMAKE_OUTPUT_CACHE" >> ~/.bashrc; \
  echo "export SNAKEMAKE_CONDA_PREFIX=$SNAKEMAKE_CONDA_PREFIX" >> ~/.bashrc

ENV PATH=/opt/conda/envs/snakemake/bin:$PATH
#################### /snakemake ####################


#################### copier ####################
FROM snakemake AS copier

RUN set -eux; \
  apt-get update; \
  apt-get install -y --no-install-recommends --no-install-suggests libc6-dev gcc; \
  pip install --no-cache-dir copier==5.1.0; \
  apt-get purge -y libc6-dev gcc; \
  apt autoremove -y; \
  apt-get clean all; \
  rm -rf /var/lib/apt/lists/* /usr/share/doc/* /usr/share/man/* /opt/conda/man/*
#################### /copier ####################


#################### workflow ####################
FROM copier AS workflow

WORKDIR /usr/src/code

COPY . /usr/src/code

VOLUME ["$SNAKEMAKE_OUTPUT_CACHE", \
        "$SNAKEMAKE_CONDA_PREFIX"]

CMD ["snakemake", "--cores=all", "--use-conda", "-p"]
#################### workflow ####################
