FROM ubuntu:20.04

# modified from https://github.com/RTIInternational/code_docker_lib/blob/master/metaxcan/Dockerfile

# Don't print "debconf: unable to initialize frontend: Dialog" messages
ARG DEBIAN_FRONTEND=noninteractive

# Install basic stuff and R
RUN apt-get update && apt-get install -y \
    locales \
    git \
    vim-tiny \
    less \
    wget \
    r-base \
    fonts-texgyre \
    texinfo

RUN locale-gen en_US.utf8 \
    && /usr/sbin/update-locale LANG=en_US.UTF-8
ENV LANG=en_US.UTF-8

RUN echo 'options(\n\
  repos = c(CRAN = "https://cloud.r-project.org/"),\n\
  download.file.method = "libcurl",\n\
  # Detect number of physical cores\n\
  Ncpus = parallel::detectCores(logical=FALSE)\n\
)' >> /etc/R/Rprofile.site


# Install R packages
RUN R -e 'install.packages(c("ggplot2", "dplyr"), repos="http://cran.us.r-project.org")'

# =====================================================================
# Install Miniconda for python dependencies
# =====================================================================
# Install basic stuff to be able to install Miniconda
RUN apt-get update && apt-get install -y bzip2 \
    ca-certificates \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender1 \
    mercurial \
    subversion

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# Set path to conda
ENV PATH /opt/conda/bin:$PATH

# Updating Anaconda packages
RUN conda update conda
RUN conda update --all

# =====================================================================
# Install Metaxcan and download sample data
# =====================================================================
RUN git clone https://github.com/hakyimlab/MetaXcan

# copy sample data
RUN mkdir /ref
COPY sample_data.tar.gz  /ref/sample_data.tar.gz
RUN tar xvf /ref/sample_data.tar.gz -C /ref && \
    chmod -R u+x /ref

# =====================================================================
# Standard docker image stuff
# =====================================================================
# Add metaxcan to /opt folder
ADD run.sh /
RUN mkdir /opt/code_docker_lib
RUN mv MetaXcan /opt/code_docker_lib
RUN mv run.sh /opt/code_docker_lib

# Change permissions to make things exectuable
RUN chmod u+x /opt/code_docker_lib/run.sh
RUN chmod -R u+x /opt/code_docker_lib/MetaXcan

# =====================================================================
# Install required conda environment
# =====================================================================
RUN conda env create -f /opt/code_docker_lib/MetaXcan/software/conda_env.yaml

# Set wrkdir
RUN mkdir -p /data
WORKDIR /data

