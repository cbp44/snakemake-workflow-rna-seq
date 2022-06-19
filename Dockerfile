FROM docker.io/snakemake/snakemake:v7.8.2 AS base

#################### base ####################
# Important to set this for using sort commands
ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# Add useful aliases for use inside the container
RUN echo "alias ls='/bin/ls -lha --color=auto --group-directories-first'" >> /etc/bash.bashrc

# Where the workflow code will be stored
ARG WORKFLOW_PATH
ENV WORKFLOW_PATH=${WORKFLOW_PATH:-/mnt/workflow}

# Location of snakemake output cache
ENV SNAKEMAKE_OUTPUT_CACHE=/mnt/snakemake-cache

# Where snakemake conda environments get created
ENV SNAKEMAKE_CONDA_PREFIX=/conda-envs
#################### /base ####################


#################### workflow ####################
FROM base AS workflow

VOLUME ["${SNAKEMAKE_OUTPUT_CACHE}", "${SNAKEMAKE_CONDA_PREFIX}"]

# Define volume for the workflow results
VOLUME ["${WORKFLOW_PATH}/resources","${WORKFLOW_PATH}/results"]

COPY ./ "${WORKFLOW_PATH}"

WORKDIR "${WORKFLOW_PATH}"

# Default command runs 3 rules at a time, downloads at most 2 at the same time
#   and uses conda for software environments.
CMD ["snakemake","--cores=all", "--printshellcmds", "--use-conda"]
#################### /workflow ###################



#################### copier ####################
# Install copier into the conda env named copier
FROM base AS copier

ENV PIP_NO_CACHE_DIR=1

RUN mamba create -y -q -n copier pip \
      && mamba clean --all -y

ENV PATH /opt/conda/envs/copier/bin:${PATH}

RUN install_packages build-essential \
      # && /opt/conda/envs/copier/bin/pip install copier==6.1.0 jinja2-time==0.2.0 \
      && /opt/conda/envs/copier/bin/pip install --no-cache-dir copier==6.1.0 copier-templates-extensions==0.2.0 \
      && apt autoremove -y build-essential

RUN echo "source activate copier" >> ~/.bashrc \
      && git config --global user.email "copier@dockercontainer.com" \
      && git config --global user.name "copier"

COPY --from=workflow "${WORKFLOW_PATH}/" ${WORKFLOW_PATH}

ENTRYPOINT ["/opt/conda/envs/copier/bin/copier"]
CMD ["--help"]
#################### /copier ####################