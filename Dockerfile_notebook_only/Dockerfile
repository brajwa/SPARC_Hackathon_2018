# this build contains just the Jupyter notebook
# Dockerfile for the Spat-An-App Project Team

# use Ubuntu 18.04 base image
FROM ubuntu:18.04

MAINTAINER cchorn "chorn@pitt.edu"

# run as root user
USER root

ENV DEBIAN_FRONTEND noninteractive

# ----------------------------------------------------------------
# set up oSparc env variables
ENV SIMCORE_NODE_UUID="-1" \
    SIMCORE_USER_ID="-1" \
    SIMCORE_NODE_BASEPATH="-1" \
    STORAGE_ENDPOINT="=1" \
    S3_ENDPOINT="=1" \
    S3_ACCESS_KEY="-1" \
    S3_SECRET_KEY="-1" \
    S3_BUCKET_NAME="-1" \
    POSTGRES_ENDPOINT="-1" \
    POSTGRES_USER="-1" \
    POSTGRES_PASSWORD="-1" \
    POSTGRES_DB="-1"

# Install OS packages
RUN apt-get update && apt-get install -y \
    apt-utils \
    wget \
    bzip2 \
    ca-certificates \
    sudo \
    locales \
    fonts-liberation \
    git \
    nano \
    libzmq3-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    virtualenv \
    curl 

RUN echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && \
    locale-gen

# Configure environment
ENV SHELL=/bin/bash \
    Notebook_USER=datawrangler \
    LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    LANGUAGE=en_US.UTF-8 \
    PYTHONPATH=/route

# Copy example notebooks to tmp directory
COPY Jupyter_Notebooks/ /home/
COPY input_retriever.py /route/

# Install Python 3 packages
####################################################

# Install python3 and pip package manager
RUN apt-get update \
  && apt-get install -y python3-pip python3-dev \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 install --upgrade pip 

# Install python3 packages ... removed rpy2 because of installation error
RUN pip3 install jupyter opencv-contrib-python pandas numpy matplotlib ipywidgets scipy notedown statsmodels seaborn blackfynn scikit-image

# Install simcore packagek
RUN pip3 install git+https://github.com/ITISFoundation/osparc-simcore.git#egg=version_subpkg&subdirectory=tree/master/packages/simcore-sdk
RUN pip3 install git+https://github.com/ITISFoundation/osparc-simcore.git#egg=version_subpkg&subdirectory=services/storage/client-sdk/python

# Install nodejss
# https://github.com/nodesource/distributions
RUN curl -sL https://deb.nodesource.com/setup_10.x | sudo -E bash -
RUN sudo apt-get install -y nodejs 

# Install yarn
RUN npm install -g yarn

# Install R and packages
####################################################

# Install R
RUN apt-get install -y r-base-dev \
    && apt-get clean \
	&& apt-get remove \
	&& rm -rf /var/lib/apt/lists/*

# Set default R CRAN repo
RUN echo 'options("repos"="http://cran.rstudio.com")' >> /usr/lib/R/etc/Rprofile.site

# Install R Packages and kernel for Jupyter notebook (not sure which R packages at this point)
RUN Rscript -e "install.packages(c('devtools', 'magrittr', 'dplyr', 'spatstat', 'png', 'httr', 'lctools', 'dixon', 'tibble', 'multcomp', 'emmeans', 'readxl', 'ggfortify'))"
RUN Rscript -e "devtools::install_github('IRkernel/IRkernel')"
RUN Rscript -e "IRkernel::installspec()"

# Configure Jupyter notebook
####################################################

RUN jupyter notebook --generate-config && \
    ipython profile create
# TextFileContentsManager is needed to jupytext
RUN echo "c.NotebookApp.open_browser = False" >>/root/.jupyter/jupyter_notebook_config.py && \
    echo "c.InteractiveShellApp.matplotlib = 'inline'" >>/root/.ipython/profile_default/ipython_config.py
        
label io.simcore.key = '{"key": "simcore/services/dynamic/spat-an-app"}' \
      && io.simcore.version = '{"version": "1.0.0"}' \
      && io.simcore.type = '{"type": "dynamic"}' \
      && io.simcore.name = '{"name": "spat-an-app"}' \
      && io.simcore.description = '{"description": "Jupyter lab"}' \
      && io.simcore.authors = '{"authors": [{"name": "cchorn", "email": "chorn@pitt.edu", "affiliation": "Univ Pittsburgh"}]}' \
      && io.simcore.contact = '{"contact": "chorn@pitt.edu"}' \
      && io.simcore.inputs = '{"inputs": {"input_1":{"label": "input 1", "displayOrder":0, "description": "Input 1", "type": "data:*/*"}, "input_2":{"label": "input 2", "displayOrder":1, "description": "Input 2", "type": "data:*/*"}}}' \
      && io.simcore.outputs = '{"outputs": {"output_1":{"label": "output 1", "displayOrder":0, "description": "Output 1", "type": "data:*/*"}, "output_2":{"label": "output 2", "displayOrder":1, "description": "Output 2", "type": "data:*/*"}}}' \
      && simcore.service.settings = '[{"name": "resources", "type": "Resources", "value": {"mem_limit":17179869184, "cpu_limit": 4000000000}}, {"name": "ports", "type": "int", "value": 8888}, {"name": "constraints", "type": "string", "value": ["node.platform.os == linux"]}]' \
      && org.label-schema.schema-version = "1.0" \
      && org.label-schema.build-date = "${BUILD_DATE}" \
      && org.label-schema.vcs-url = "https://github.com/ITISFoundation/osparc-simcore" \
      && org.label-schema.vcs-ref = "${VCS_REF}"

# Run the Jupyter lab .. comment the first command because it is only for the notebook
CMD jupyter notebook --allow-root --ip 0.0.0.0 --no-browser --NotebookApp.nbserver_extensions="{'input_retriever':True}"

# set directories and ports
####################################################

WORKDIR /home
RUN pwd

# Expose port to host
EXPOSE 8888
