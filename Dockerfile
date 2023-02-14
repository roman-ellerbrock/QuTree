FROM debian:stable-slim AS compile-image
RUN apt-get update
RUN apt-get install -y --no-install-recommends \
        build-essential gcc gfortran cmake \
	libboost-all-dev \
	libyaml-cpp-dev \
	liblapacke-dev

WORKDIR /app-build
COPY . ./

RUN cmake -S . -B build -DCMAKE_PREFIX_PATH=/opt/qutree -DCMAKE_INSTALL_PREFIX=/opt/qutree -DGITHUB_ACTIONS=true
RUN cmake --build build -j$(nproc) --target QuTree Hamiltonians mctdh
RUN cmake --install build
RUN cmake --install build/contrib



FROM debian:stable-slim AS build-image
COPY --from=compile-image /opt/qutree /opt/qutree

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
	libyaml-cpp-dev \
        liblapacke-dev && \
    apt-get clean && apt-get autoremove

ENV PATH="/opt/qutree/bin:$PATH"
ENV LD_LIBRARY_PATH="/opt/qutree/lib:$LD_LIBRARY_PATH"
WORKDIR /app
