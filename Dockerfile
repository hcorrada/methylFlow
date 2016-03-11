FROM ubuntu:latest

RUN apt-get update && apt-get install -y \
  gcc \
  g++ \
  cmake

COPY . /usr/src/methylFlow
WORKDIR /usr/src/methylFlow

RUN rm -rf build && \
  mkdir build && \
  cd build && \
  cmake .. && \
  make && \
  make test
