FROM ubuntu:18.04
RUN apt-get update -y && apt-get install -y autoconf build-essential libboost-all-dev zlib1g-dev libbz2-dev libcurl4-openssl-dev liblzma-dev samtools cmake python3 python3-pip libjpeg-dev
RUN pip3 install --upgrade pip setuptools wheel
RUN pip3 install matplotlib
RUN apt-get install -y python curl bzip2 && apt-get autoclean && rm -rf /var/lib/apt/lists/*
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf -
ENV PATH="${PATH}:/minimap2-2.17_x64-linux"
COPY . /workdir
WORKDIR /workdir
RUN mkdir build && cd build && cmake .. && make -j 4
ENV PATH="${PATH}:/workdir/build"