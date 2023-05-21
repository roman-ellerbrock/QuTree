# Set the base image
FROM ubuntu:latest

ARG DEBIAN_FRONTEND=noninteractive

# Update and install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    wget \
    unzip

RUN apt-get install -y \
    # yaml 
    libyaml-cpp-dev \
    # for googletests
    git \
    # for debugging
    vim

# Download LibTorch (check the version and cuda compatibility)
# Potential conflicts with yaml-cpp: https://github.com/pytorch/pytorch/issues/19353
# current linux repo for CPU & >=C++11, select at https://pytorch.org/get-started/locally/
# ...and unzip the correct file! Note: %2B in download link translates to "+" locally
RUN wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcpu.zip
RUN unzip libtorch-cxx11-abi-shared-with-deps-2.0.1+cpu.zip

# Set environment variables for LibTorch
ENV TORCH_HOME=/libtorch
ENV LD_LIBRARY_PATH=/libtorch/lib:$LD_LIBRARY_PATH
ENV PATH=$PATH:/libtorch/bin

# Copy your project into the Docker image
COPY . /qutree

# Go into the project directory
WORKDIR /qutree

# Build your project
RUN mkdir build

WORKDIR /tetrachem/build

RUN cmake -DCMAKE_PREFIX_PATH=/libtorch ..

RUN make -j 2

WORKDIR /qutree/build/tests

#CMD ["./qutree/build/tests/tests"]