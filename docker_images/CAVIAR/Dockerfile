# Note that this container will only work on x86_64 systems due to qemu dependency
FROM amd64/ubuntu:jammy
LABEL maintainer="Wendy Wong <wendy.wong@gmail.com>"    \
      software="caviar" 

RUN apt-get update && apt install -y libgsl-dev liblapack-dev libblas-dev git make g++-9  && apt clean 
RUN ln -s /usr/bin/g++-9 /usr/bin/g++
RUN git clone https://github.com/fhormoz/caviar.git && cd caviar/CAVIAR-C++/ && make && cp *CAVIAR /usr/local/bin/
ENV PATH="/usr/local/bin/:${PATH}"