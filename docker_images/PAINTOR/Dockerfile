FROM amd64/ubuntu:jammy
LABEL maintainer="Wendy Wong <wendy.wong@gmail.com>"    \
      software="PAINTOR_V3.0" 

RUN apt-get update && apt install -y wget make g++ && apt clean 
RUN wget https://github.com/gkichaev/PAINTOR_V3.0/archive/refs/tags/3.0.tar.gz && tar -zxvf 3.0.tar.gz && cd PAINTOR_V3.0-3.0 && bash install.sh && cp PAINTOR /usr/local/bin/