FROM ubuntu:20.04

RUN apt-get update && apt-get install -y software-properties-common gcc && \
    add-apt-repository -y ppa:deadsnakes/ppa

RUN apt-get update && apt-get install -y python3.8 python3-pip

RUN apt-get install -y build-essential g++

RUN apt-get install -y python3-distutils

RUN pip install 'scanpy[leiden]'

RUN pip install numba

RUN apt-get install -y wget

COPY scanpy_processing.py /usr/local/bin/

RUN chmod a+rx /usr/local/bin/scanpy*.py

ENV PATH="${PATH}:/usr/local/bin"