FROM ubuntu:16.04

LABEL org.opencontainers.image.source="https://github.com/guillaume-gricourt/HmnTrimmer"

RUN apt-get update
RUN apt-get install -y \
    build-essential \
    python3 \
    wget \
    yasm \
    zlib1g-dev

# libssl-dev libffi-dev \
RUN wget https://bootstrap.pypa.io/pip/3.5/get-pip.py
RUN python3 get-pip.py --ignore-installed
RUN pip3 install django matplotlib seaborn packaging

WORKDIR /opt
RUN mkdir HmnTrimmer
COPY . /opt/HmnTrimmer/

RUN cd /opt/HmnTrimmer/lib/igzip-042/igzip && \
    make slib0c
RUN cd /opt/HmnTrimmer && \
    make
RUN ln /opt/HmnTrimmer/HmnTrimmer /usr/local/bin/
ENV LD_LIBRARY_PATH=/opt/HmnTrimmer/lib/igzip-042/igzip

ENTRYPOINT ["HmnTrimmer"]
