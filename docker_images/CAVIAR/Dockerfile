FROM ubuntu:focal
LABEL maintainer="Wendy Wong <wendy.wong@gmail.com>"    \
      software="caviar" 

RUN apt-get update && apt install -y libgsl-dev liblapack-dev libblas-dev git make g++
RUN git clone https://github.com/fhormoz/caviar.git && cd caviar/CAVIAR-C++/ && make && cp *CAVIAR /usr/local/bin/
ENV PATH="/usr/local/bin/:${PATH}"