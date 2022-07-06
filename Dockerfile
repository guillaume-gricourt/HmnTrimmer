FROM ubuntu:16.04

LABEL org.opencontainers.image.source="https://github.com/guillaume-gricourt/HmnTrimmer"

RUN apt-get update
RUN apt-get install -y \
    build-essential \
    python3 \
    wget \
    yasm \
    zlib1g-dev

WORKDIR /opt
# Install python
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    chmod +x Miniconda3-latest-Linux-x86_64.sh && \
    ./Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda3 && \
    rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH="/opt/miniconda3/bin/:$PATH"
# Install deps
RUN /opt/miniconda3/bin/pip3 install django matplotlib seaborn packaging

# Install HmnTrimmer
RUN mkdir HmnTrimmer
COPY . /opt/HmnTrimmer/

RUN cd /opt/HmnTrimmer/lib/igzip-042/igzip && \
    make slib0c
RUN cd /opt/HmnTrimmer && \
    make
RUN ln /opt/HmnTrimmer/HmnTrimmer /usr/local/bin/
ENV LD_LIBRARY_PATH=/opt/HmnTrimmer/lib/igzip-042/igzip

ENTRYPOINT ["HmnTrimmer"]
