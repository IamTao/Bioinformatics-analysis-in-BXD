## Based on
FROM taptao/rstudio

MAINTAINER itamtao@gmail.com

## install some necessary tools
RUN apt-get -y update && apt-get install -y curl
RUN apt-get update && apt-get install -y zsh
RUN apt-get install -y git

RUN apt-get install -y --no-install-recommends \
    python \
    python-dev

RUN curl -O https://bootstrap.pypa.io/get-pip.py && \
    python get-pip.py && \
    rm get-pip.py

## install the required R packages

## install the required python packages
RUN pip install --no-cache-dir numpy
RUN pip --no-cache-dir install pandas
