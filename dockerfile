FROM ubuntu:16.04

LABEL application="HmnTrimmer"
LABEL maintainer="guillaume.gricourt@aphp.fr"

RUN apt-get update && \
    apt-get install -y build-essential && \
    apt-get install -y zlib1g-dev

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
