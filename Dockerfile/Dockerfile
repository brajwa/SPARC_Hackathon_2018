# Dockerfile for the Spat-An-App Project Team

# use Ubuntu 18.04 base image
FROM ubuntu:18.04

MAINTAINER cchorn "chorn@pitt.edu"

# run as root user
USER root

ENV DEBIAN_FRONTEND noninteractive

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
    LANGUAGE=en_US.UTF-8

# Copy example notebooks to tmp directory
COPY Jupyter_Notebooks/ /home/

# Install Python 3 packages
####################################################

# Install python3 and pip package manager
RUN apt-get update \
  && apt-get install -y python3-pip python3-dev \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 install --upgrade pip 

# Install python3 packages 
RUN pip3 install jupyter jupyterlab pandas numpy matplotlib ipywidgets scipy notedown statsmodels seaborn blackfynn scikit-image rpy2==2.8.5

# Install Jupyter lab extensions
####################################################

# Install nodejs for Jupyter lab extensions
# https://github.com/nodesource/distributions
RUN curl -sL https://deb.nodesource.com/setup_10.x | sudo -E bash -
RUN sudo apt-get install -y nodejs 

# Install yarn
RUN npm install -g yarn

# Install Jupyterlab extensions
RUN jupyter labextension install @jupyterlab/toc @jupyter-widgets/jupyterlab-manager

# Install R and packages
####################################################

# Install R
RUN apt-get install r-base-dev \
	&& apt-get clean \
	&& apt-get remove \
	&& rm -rf /var/lib/apt/lists/*

# Set default R CRAN repo
RUN echo 'options("repos"="http://cran.rstudio.com")' >> /usr/lib/R/etc/Rprofile.site

# Install R Packages and kernel for Jupyter notebook (not sure which R packages at this point)
RUN Rscript -e "install.packages(c('devtools'))"
RUN Rscript -e "devtools::install_github('IRkernel/IRkernel')"
RUN Rscript -e "IRkernel::installspec()"

# Configure Jupyter notebook
####################################################

RUN jupyter notebook --generate-config && \
    ipython profile create
# TextFileContentsManager is needed to jupytext
RUN echo "c.NotebookApp.open_browser = False" >>/root/.jupyter/jupyter_notebook_config.py && \
    echo "c.InteractiveShellApp.matplotlib = 'inline'" >>/root/.ipython/profile_default/ipython_config.py && \
    echo "c.NotebookApp.contents_manager_class = 'jupytext.TextFileContentsManager'" >>/root/.jupyter/jupyter_notebook_config.py

# Run the Jupyter lab .. comment the first command because it is only for the notebook
CMD jupyter lab --allow-root --ip 0.0.0.0 --no-browser

# set directories and ports
####################################################

WORKDIR /home
RUN pwd

# Expose port to host
EXPOSE 8888
